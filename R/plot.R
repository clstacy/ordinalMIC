utils::globalVariables(c(
  # results tables
  "MIC","CI_Lower","CI_Upper","Delta_MIC","Ratio_MIC",
  "Group1","Group2","DDlog2MIC","x_simplified","facet",
  # observed MIC assembly / plotting
  "._conc_raw_","._rep_","._key_","._bar_pos_",
  "x_val","color_val","y","y_obs","y_cens","cens",
  "x_center","color_off","rep_offset","x_pos","gkey",
  # misc
  ".data","pretty_label","rep","censored","max_tested","visible"
))


#' autoplot method for mic_solve objects
#'
#' Quick ggplot2 visualisations of the main outputs from
#' [mic_solve()].  Three panels are supported:
#' * **"mic"**   - forest plot of group-wise MIC estimates with
#'   asymmetric CIs.
#' * **"delta"** - forest plot of deltaMIC pairwise differences.
#' * **"ratio"** - forest plot of MIC ratios (log scale).
#'
#' @param x        Variable for x axis plotting
#' @param object   An object returned by [mic_solve()].
#' @param type     One of `"mic"`, `"delta"`, or `"ratio"`, `"DoD_delta"`, or `"DoD_ratio"`.
#' @param color_by Optional column name used to color and dodge replicate
#'   points. Default: first column in `newdata`.
#' @param dot_size Size of the dots in the dotplot. Default: `0.5`.
#' @param ...      Additional arguments passed to [ggplot2::ggplot()].
#'
#' @return A **ggplot** object.
#' @author Carson Stacy
#' @examples
#' if (requireNamespace("ordinal", quietly = TRUE) &&
#'     requireNamespace("ggplot2", quietly = TRUE)) {
#'   df <- data.frame(score  = ordered(sample(0:4, 120, TRUE)),
#'                    conc   = runif(120, 0, 4),
#'                    strain = factor(sample(c("A","B"), 120, TRUE)))
#'   fit <- ordinal::clm(score ~ strain * log1p(conc), data = df)
#'   res <- mic_solve(fit, expand.grid(strain = levels(df$strain)),
#'                    conc_name = "conc")
#'   ggplot2::autoplot(res, type = "mic")
#' }
#'
#' @name autoplot.mic_solve
#' @method autoplot mic_solve
#' @exportS3Method ggplot2::autoplot mic_solve
#' @aliases plot.mic_solve
#' @importFrom ggplot2 autoplot
#' @importFrom ggplot2 ggplot aes geom_pointrange geom_errorbar theme_minimal labs coord_flip theme element_text scale_y_log10
#' @importFrom stats ave relevel setNames formula terms
#'
autoplot.mic_solve <- function(object,
                               type = c("mic", "delta", "ratio",
                                        "DoD_delta", "DoD_ratio"),
                               x = NULL, color_by = NULL,# shape_by = NULL,
                               dot_size=0.5,
                               ...) {
  type <- match.arg(type)
  stopifnot(inherits(object, "mic_solve"))
  pd <- ggplot2::position_dodge(width = 0.5)

  if (type == "mic") {
    ## ---- inputs & aesthetics ----------------------------------------------
    tbl <- object$mic_estimates
    group_vars <- names(object$newdata)
    if (is.null(x)) x <- group_vars[1]
    if (is.null(color_by)) color_by <- group_vars[min(2, length(group_vars))]

    # keep ordering consistent with original data
    x_lvls <- if (x %in% names(object$data)) levels(object$data[[x]]) else unique(tbl[[x]])
    if (is.null(x_lvls)) x_lvls <- unique(tbl[[x]])
    col_lvls <- if (color_by %in% names(object$data)) levels(object$data[[color_by]]) else unique(tbl[[color_by]])
    if (is.null(col_lvls)) col_lvls <- "(1)"

    # base centers per x group
    x_index <- setNames(seq_along(x_lvls), x_lvls)

    # color-slice dodge across the x group
    dodge_w_colors <- 0.6
    M <- length(col_lvls)
    color_offsets <- setNames(
      if (M == 1) 0 else ((seq_len(M) - (M + 1)/2) / max(1, (M/2))) * (dodge_w_colors/2),
      col_lvls
    )
    color_band_w <- if (M == 1) dodge_w_colors else (dodge_w_colors / M)

    # add continuous x pos for bars
    tbl[[x]] <- factor(tbl[[x]], levels = x_lvls)
    if (color_by %in% names(tbl)) tbl[[color_by]] <- factor(tbl[[color_by]], levels = col_lvls)
    tbl$._bar_pos_ <- unname(x_index[as.character(tbl[[x]])]) +
      if (color_by %in% names(tbl)) unname(color_offsets[as.character(tbl[[color_by]])]) else 0

    ## ---- observed MIC per (group x replicate series) -----------------------
    d <- as.data.frame(object$data)
    score_col <- names(object$model$model)[1]

    # original-scale concentration for plotting (avoid double inverse)
    has_trans_col <- (!identical(object$transform_fun, identity)) &&
      (object$trans_name %in% names(d))

    if (identical(object$transform_fun, identity)) {
      d$._conc_raw_ <- suppressWarnings(as.numeric(d[[object$conc_name]]))
    } else if (has_trans_col) {
      d$._conc_raw_ <- object$inv_transform_fun(
        suppressWarnings(as.numeric(d[[object$trans_name]]))
      )
    } else {
      # inline transform; no transformed column stored -> use raw as-is
      d$._conc_raw_ <- suppressWarnings(as.numeric(d[[object$conc_name]]))
    }

    # ensure a replicate column (series across conc)
    rep_col <- if ("rep" %in% names(d)) "rep" else NULL
    if (is.null(rep_col)) {
      conc_col_used <- if (has_trans_col) object$trans_name else object$conc_name
      within_gc <- interaction(d[c(group_vars, conc_col_used)], drop = TRUE)
      # Correctly generate replicate numbers within each group
      d$._rep_ <- ave(rep(1, nrow(d)), within_gc, FUN = seq_along)
      rep_col <- "._rep_"
    }

    # split series by (x, color_by, replicate)
    split_keys <- unique(na.omit(c(x, if (color_by %in% names(d)) color_by, rep_col)))
    d$._key_ <- interaction(d[split_keys], drop = TRUE)
    chunks <- split(d, d$._key_, drop = TRUE)

    # global baseline from the model frame (same scale as sc below)
    global_base <- suppressWarnings(
      min(as.numeric(as.character(object$model$model[[1]])), na.rm = TRUE)
    )

    .as_num <- function(v) {
      out <- suppressWarnings(as.numeric(as.character(v)))
      if (all(is.na(out))) as.numeric(v) else out
    }

    obs_rows <- lapply(chunks, function(gdf) {
      x_val     <- as.character(gdf[[x]][1])
      color_val <- if (color_by %in% names(gdf)) as.character(gdf[[color_by]][1]) else "(1)"

      # order by raw concentration
      conc <- suppressWarnings(as.numeric(gdf$._conc_raw_))
      ord  <- order(conc)
      conc <- conc[ord]
      sc   <- .as_num(gdf[[score_col]][ord])

      # collapse duplicate rows at the SAME concentration within this series:
      # mark 'growth at conc' if ANY row at that conc has sc > 0 (global base)
      grew_by_conc <- tapply(sc > global_base, conc, function(v) any(v, na.rm = TRUE))
      cvals <- as.numeric(names(grew_by_conc))
      grew  <- as.logical(grew_by_conc)
      o2    <- order(cvals)
      cvals <- cvals[o2]
      grew  <- grew[o2]

      n <- length(cvals)
      if (n == 0L) {
        return(data.frame(
          x_val     = x_val,
          color_val = color_val,
          y_obs     = NA_real_,
          y_cens    = NA_real_,
          stringsAsFactors = FALSE
        ))
      }

      if (any(grew)) {
        last_g <- max(which(grew))
        if (last_g < n) {
          mic_val <- cvals[last_g + 1L]   # earliest conc AFTER last growth
          y_obs   <- mic_val
          y_cens  <- NA_real_
        } else {
          # growth persists to the highest tested concentration -> censored
          y_obs  <- NA_real_
          y_cens <- max(cvals, na.rm = TRUE)
        }
      } else {
        # no growth anywhere -> MIC is the lowest tested concentration
        y_obs  <- cvals[1L]
        y_cens <- NA_real_
      }

      data.frame(
        x_val     = x_val,
        color_val = color_val,
        y_obs     = y_obs,
        y_cens    = y_cens,
        stringsAsFactors = FALSE
      )
    })


    obs_tbl <- do.call(rbind, obs_rows)
    if (!nrow(obs_tbl)) {
      # draw only bars if no points
      bar_width <- color_band_w * 0.85
      p <- ggplot2::ggplot() +
        ggplot2::geom_col(
          data = tbl,
          mapping = ggplot2::aes(x = ._bar_pos_, y = MIC,
                                 fill = if (color_by %in% names(tbl)) .data[[color_by]] else NULL,
                                 color = if (color_by %in% names(tbl)) .data[[color_by]] else NULL),
          width = bar_width, alpha = 0.5, color = NA
        ) +
        ggplot2::geom_errorbar(
          data = tbl,
          mapping = ggplot2::aes(x = ._bar_pos_, ymin = CI_Lower, ymax = CI_Upper,
                                 color = if (color_by %in% names(tbl)) .data[[color_by]] else NULL),
          width = bar_width * 0.25
        ) +
        ggplot2::coord_flip() +
        ggplot2::scale_x_continuous(breaks = seq_along(x_lvls), labels = x_lvls,
                                    expand = ggplot2::expansion(mult = c(0.05, 0.1))) +
        ggplot2::theme_minimal() +
        ggplot2::labs(x = "Group", y = "MIC", fill = color_by, color = color_by,
                      title = "Modelled MIC* with observed MICs") +
        ggplot2::theme(legend.position = "bottom")
      return(p)
    }

    obs_tbl$x_val     <- factor(obs_tbl$x_val, levels = x_lvls)
    obs_tbl$color_val <- factor(obs_tbl$color_val, levels = col_lvls)

    ## ---- “dotplot-style” replicate offsets per (x,color,y) ------------------
    slice_w <- color_band_w * 0.3#7

    # precompute base centers
    obs_tbl$x_center  <- unname(x_index[as.character(obs_tbl$x_val)])
    obs_tbl$color_off <- unname(color_offsets[as.character(obs_tbl$color_val)])

    # combine observed/censored into one table with y + cens flag
    # make_pts <- function(df, ycol, cens_flag) {
    #   if (!nrow(df)) {
    #     return(data.frame(
    #       x_val=character(0), color_val=character(0),
    #       y=numeric(0), cens=logical(0),
    #       stringsAsFactors = FALSE
    #     ))
    #   }
    #   data.frame(
    #     x_val=df$x_val, color_val=df$color_val,
    #     y=df[[ycol]], cens=rep(cens_flag, nrow(df)),
    #     stringsAsFactors = FALSE
    #   )
    # }
#
#     obs_ok  <- subset(obs_tbl, is.finite(y_obs))
#     cens_ok <- subset(obs_tbl, is.finite(y_cens))
#
#     pts_all <- rbind(
#       make_pts(obs_ok,  "y_obs",  FALSE),
#       make_pts(cens_ok, "y_cens", TRUE)
#     )
#
#     pts_all$gkey <- interaction(pts_all$x_val, pts_all$color_val, pts_all$y, drop = TRUE)

    # --- safe builders for observed/censored point tables ---
    make_pts <- function(df, ycol, cens_flag) {
      if (!nrow(df)) {
        return(data.frame(
          x_val=character(0), color_val=character(0),
          y=numeric(0), cens=logical(0),
          stringsAsFactors = FALSE
        ))
      }
      data.frame(
        x_val=df$x_val, color_val=df$color_val,
        y=df[[ycol]], cens=rep(cens_flag, nrow(df)),
        stringsAsFactors = FALSE
      )
    }

    obs_ok  <- subset(obs_tbl, is.finite(y_obs))
    cens_ok <- subset(obs_tbl, is.finite(y_cens))

    pts_all <- rbind(
      make_pts(obs_ok,  "y_obs",  FALSE),
      make_pts(cens_ok, "y_cens", TRUE)
    )

    # if no points, bail out to bars-only plot (same as your earlier early-return)
    if (!nrow(pts_all)) {
      bar_width <- color_band_w * 0.85
      p <- ggplot2::ggplot() +
        ggplot2::geom_col(
          data = tbl,
          mapping = ggplot2::aes(x = ._bar_pos_, y = MIC,
                                 fill  = if (color_by %in% names(tbl)) .data[[color_by]] else NULL,
                                 color = if (color_by %in% names(tbl)) .data[[color_by]] else NULL),
          width = bar_width, alpha = 0.5, color = NA
        ) +
        ggplot2::geom_errorbar(
          data = tbl,
          mapping = ggplot2::aes(x = ._bar_pos_, ymin = CI_Lower, ymax = CI_Upper,
                                 color = if (color_by %in% names(tbl)) .data[[color_by]] else NULL),
          width = bar_width * 0.25
        ) +
        ggplot2::coord_flip() +
        ggplot2::scale_x_continuous(breaks = seq_along(x_lvls), labels = x_lvls,
                                    expand = ggplot2::expansion(mult = c(0.05, 0.1))) +
        ggplot2::theme_minimal() +
        ggplot2::labs(x = "Group", y = "MIC", fill = color_by, color = color_by,
                      title = "Modelled MIC* with observed MICs") +
        ggplot2::theme(legend.position = "bottom")
      return(p)
    }

    # --- reattach position columns for dot placement ---
    # (use the same lookups you used to build obs_tbl)
    pts_all$x_val    <- factor(pts_all$x_val, levels = x_lvls)
    pts_all$color_val<- factor(pts_all$color_val, levels = col_lvls)

    pts_all$x_center  <- unname(x_index[as.character(pts_all$x_val)])
    pts_all$color_off <- if (length(color_offsets)) {
      co <- color_offsets[as.character(pts_all$color_val)]
      unname(ifelse(is.na(co), 0, co))
    } else 0

    # cluster key for dotplot-like offsetting
    pts_all$gkey <- interaction(pts_all$x_val, pts_all$color_val, pts_all$y, drop = TRUE)

    # helpers to generate dotplot-like offsets (unchanged)
    make_offsets <- function(L, width, prefer_center) {
      if (L == 1) return(0)
      if (prefer_center && (L %% 2 == 1)) {
        k <- (L - 1L) / 2L
        raw <- c(-seq_len(k), 0, seq_len(k))
        return((raw / max(abs(raw))) * (slice_w / 2))
      } else {
        steps <- seq(0.5, by = 1, length.out = ceiling(L / 2))
        raw   <- sort(c(-steps, steps))[seq_len(L)]
        return((raw / max(abs(raw))) * (slice_w / 2))
      }
    }

    pts_all$rep_offset <- NA_real_
    cluster_index <- split(seq_len(nrow(pts_all)), pts_all$gkey)

    for (idx in cluster_index) {
      L <- length(idx)
      has_obs_idx <- which(!pts_all$cens[idx])
      prefer_center <- length(has_obs_idx) > 0
      offs <- make_offsets(L, slice_w, prefer_center)

      if (L == 1) {
        pts_all$rep_offset[idx] <- offs
        next
      }

      if (prefer_center && (L %% 2 == 1)) {
        zero_pos <- which.min(abs(offs))
        center_i <- idx[has_obs_idx[1]]
        pts_all$rep_offset[center_i] <- offs[zero_pos]

        remaining <- setdiff(idx, center_i)
        rem_offs  <- offs[-zero_pos]
        rem_obs   <- setdiff(idx[has_obs_idx], center_i)
        rem_cens  <- setdiff(remaining, rem_obs)
        order_off <- order(abs(rem_offs), rem_offs)
        ordered_members <- c(rem_obs, rem_cens)
        pts_all$rep_offset[ordered_members] <- rem_offs[order_off][seq_along(ordered_members)]
      } else {
        order_off <- order(abs(offs), offs)
        if (length(has_obs_idx)) {
          obs_members  <- idx[has_obs_idx]
          cens_members <- setdiff(idx, obs_members)
          k <- min(length(obs_members), length(order_off))
          pts_all$rep_offset[obs_members]    <- offs[order_off][seq_len(k)]
          if (length(cens_members)) {
            pts_all$rep_offset[cens_members] <- offs[order_off][-seq_len(k)]
          }
        } else {
          pts_all$rep_offset[idx] <- offs[order_off]
        }
      }
    }

    pts_all$x_pos <- pts_all$x_center + pts_all$color_off + pts_all$rep_offset


    # helpers to generate dotplot-like offsets
    make_offsets <- function(L, width, prefer_center) {
      if (L == 1) return(0)
      if (prefer_center && (L %% 2 == 1)) {
        # odd size WITH a preferred center -> include 0
        k <- (L - 1L) / 2L
        raw <- c(-seq_len(k), 0, seq_len(k))
        return((raw / max(abs(raw))) * (width / 2))
      } else {
        # even OR no preferred center -> no 0, symmetric ±0.5, ±1.5, ...
        steps <- seq(0.5, by = 1, length.out = ceiling(L / 2))
        raw   <- sort(c(-steps, steps))[seq_len(L)]
        return((raw / max(abs(raw))) * (width / 2))
      }
    }


    pts_all$rep_offset <- NA_real_
    cluster_index <- split(seq_len(nrow(pts_all)), pts_all$gkey)

    for (idx in cluster_index) {
      L <- length(idx)
      has_obs_idx <- which(!pts_all$cens[idx])  # positions within 'idx' that are observed
      prefer_center <- length(has_obs_idx) > 0  # only “prefer center” when any observed

      offs <- make_offsets(L, slice_w, prefer_center)

      if (L == 1) {
        pts_all$rep_offset[idx] <- offs
        next
      }

      if (prefer_center && (L %% 2 == 1)) {
        # odd & with observed -> put ONE observed on 0, others by |offset|
        zero_pos <- which.min(abs(offs))              # the center position (0)
        # pick one observed to be the center
        center_i <- idx[has_obs_idx[1]]
        pts_all$rep_offset[center_i] <- offs[zero_pos]

        # remaining members -> assign offsets by increasing |offset| (skipping center)
        remaining <- setdiff(idx, center_i)
        rem_offs  <- offs[-zero_pos]
        # give observed the smallest |offset|s first, then censored
        rem_obs   <- setdiff(idx[has_obs_idx], center_i)
        rem_cens  <- setdiff(remaining, rem_obs)
        order_off <- order(abs(rem_offs), rem_offs)   # nearest to center first
        ordered_members <- c(rem_obs, rem_cens)
        pts_all$rep_offset[ordered_members] <- rem_offs[order_off][seq_along(ordered_members)]
      } else {
        # even cluster OR all-censored: no dot at 0; assign nearest-to-center slots
        # give observed the closest pair (±0.5), then fill outward
        order_off <- order(abs(offs), offs)
        if (length(has_obs_idx)) {
          obs_members  <- idx[has_obs_idx]
          cens_members <- setdiff(idx, obs_members)

          k <- length(obs_members)
          # cap in case more observed than available offsets (shouldn’t happen often)
          k <- min(k, length(order_off))

          pts_all$rep_offset[obs_members]  <- offs[order_off][seq_len(k)]
          if (length(cens_members)) {
            pts_all$rep_offset[cens_members] <- offs[order_off][-seq_len(k)]
          }
        } else {
          # all censored: just fill by closeness to center
          pts_all$rep_offset[idx] <- offs[order_off]
        }
      }
    }


    pts_all$x_pos <- pts_all$x_center + pts_all$color_off + pts_all$rep_offset


    ## ---- draw plot ----------------------------------------------------------
    bar_width <- color_band_w * 0.85
    p <- ggplot2::ggplot() +
      ggplot2::geom_col(
        data = tbl,
        mapping = ggplot2::aes(x = ._bar_pos_, y = MIC,
                               fill = if (color_by %in% names(tbl)) .data[[color_by]] else NULL,
                               color = if (color_by %in% names(tbl)) .data[[color_by]] else NULL),
        width = bar_width, alpha = 0.5, color = NA
      ) +
      ggplot2::geom_errorbar(
        data = tbl,
        mapping = ggplot2::aes(x = ._bar_pos_, ymin = CI_Lower, ymax = CI_Upper,
                               color = if (color_by %in% names(tbl)) .data[[color_by]] else NULL),
        width = bar_width * 0.25
      )

    p <- p +
      ggplot2::geom_point(
        data = pts_all,
        mapping = ggplot2::aes(
          x = x_pos, y = y,
          fill  = if (color_by %in% names(tbl)) color_val else NULL,
          color = if (color_by %in% names(tbl)) color_val else NULL,
          alpha = ifelse(cens, 0.4, 1)
        ),
        size = 2, show.legend = TRUE
      ) +
      ggplot2::scale_alpha_identity() +
      ggplot2::coord_flip() +
      ggplot2::scale_x_continuous(breaks = seq_along(x_lvls), labels = x_lvls,
                                  expand = ggplot2::expansion(mult = c(0.05, 0.1))) +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        x = "Group", y = "MIC",
        fill = color_by, color = color_by,
        title = "Modelled MIC* with observed MICs"
      ) +
      ggplot2::theme(legend.position = "bottom")

    return(p)

  } else if (type == "delta") {

    ## ----- prettier deltaMIC forest plot  ------------------------------------
    tbl <- object$delta_mic_results

    ## 1) explode Group1 / Group2 into factor columns
    fvars <- names(object$newdata)              # e.g. "strain", "treatment"
    if (is.null(color_by)) color_by <- fvars[min(2, length(fvars))]

    explode_delta <- function(col) {
      mat <- do.call(rbind, strsplit(as.character(tbl[[col]]), ":", fixed = TRUE))
      out <- as.data.frame(mat, stringsAsFactors = FALSE)
      names(out) <- fvars
      out
    }

    left  <- explode_delta("Group1")
    right <- explode_delta("Group2")
    names(left)  <- paste0(fvars, "_1")
    names(right) <- paste0(fvars, "_2")
    tbl <- cbind(tbl, left, right)

    for (i in seq_along(fvars)) {
      tbl[[fvars[i]]] <- ifelse(
        tbl[[paste0(fvars[i], "_1")]] == tbl[[paste0(fvars[i], "_2")]],
        tbl[[paste0(fvars[i], "_1")]],
        paste0(tbl[[paste0(fvars[i], "_2")]], " - ", tbl[[paste0(fvars[i], "_1")]])
      )
    }

    explode <- function(x) {
      m <- do.call(rbind, strsplit(as.character(x), ":", fixed = TRUE))
      out <- as.data.frame(m, stringsAsFactors = FALSE)
      names(out) <- fvars
      out
    }

    g1 <- explode(tbl$Group1)
    g2 <- explode(tbl$Group2)

    ## 2) facet label (put in data!)
    same <- g1 == g2
    tbl$facet <- apply(same, 1, function(r) {
      if (all(r))           "all"
      else if (!any(r))     "Effect of other"
      else                  paste0("Effect of ", paste(fvars[!r], collapse = "+"))
    })

    tbl$facet <- .relevel_if_present(tbl$facet, "Effect of other")
    tbl$facet <- factor(tbl$facet, levels = rev(levels(tbl$facet)))

    ## 3) simplified x labels
    tbl$x_simplified <- vapply(
      seq_len(nrow(tbl)),
      FUN.VALUE = character(1),
      function(i) {
        dif <- fvars[g1[i, ] != g2[i, ]]
        sam <- fvars[g1[i, ] == g2[i, ]]
        if (length(sam) == 0) sam <- NA
        if (length(dif) == 0) {
          "no-diff"
        } else if (length(dif) == 1) {
          if (is.na(sam)) {
            sprintf("%s: %s - %s", dif, g2[i, dif], g1[i, dif])
          } else {
            sprintf("%s: %s\n%s: %s - %s",
                    sam, g2[i, sam],
                    dif, g2[i, dif], g1[i, dif])
          }
        } else {
          paste(paste(g1[i, fvars], collapse = ":"),
                paste(g2[i, fvars], collapse = ":"), sep = " - ")
        }
      }
    )

    pd <- ggplot2::position_dodge(width = 0.6)

    p <- ggplot2::ggplot(
      tbl,
      ggplot2::aes(
        x = x_simplified,
        y = Delta_MIC,
        ymin = CI_Lower, ymax = CI_Upper
      )
    )

    # map color/fill via .data pronoun so it evaluates inside each layer's data
    if (!is.null(color_by) && color_by %in% names(tbl)) {
      p <- p + ggplot2::aes(color = .data[[color_by]], fill = .data[[color_by]])
    }

    p <- p +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
      ggplot2::geom_col(position = pd, width = 0.5, alpha = 0.35, color = NA) +
      ggplot2::geom_errorbar(position = pd, width = 0.25) +
      ggplot2::geom_point(position = pd, size = 2) +
      ggplot2::coord_flip() +
      ggplot2::facet_grid(rows = ggplot2::vars(facet), scales = "free_y", space = "free") +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        x = NULL, y = "delta-MIC (absolute)",
        color = color_by, fill = color_by,
        title = "Pairwise MIC differences (delta-MIC)"
      ) +
      ggplot2::theme(
        strip.text.y = ggplot2::element_text(face = "bold"),
        legend.position = if (is.null(color_by)) "none" else "bottom"
      )

    return(p)



  } else if (type == "ratio") {

    tbl <- object$ratio_mic_results
    if (is.null(tbl)) {
      stop("No pairwise MIC ratios found. Check that the model has at least two levels for each factor.")
    }

    fvars <- names(object$newdata)
    if (is.null(color_by)) color_by <- fvars[min(2, length(fvars))]

    explode_mic <- function(col) {
      mat <- do.call(rbind, strsplit(as.character(tbl[[col]]), ":", fixed = TRUE))
      out <- as.data.frame(mat, stringsAsFactors = FALSE)
      names(out) <- fvars
      out
    }

    left  <- explode_mic("Group1")
    right <- explode_mic("Group2")
    names(left)  <- paste0(fvars, "_1")
    names(right) <- paste0(fvars, "_2")
    tbl <- cbind(tbl, left, right)

    for (i in seq_along(fvars)) {
      tbl[[fvars[i]]] <- ifelse(
        tbl[[paste0(fvars[i], "_1")]] == tbl[[paste0(fvars[i], "_2")]],
        tbl[[paste0(fvars[i], "_1")]],
        paste0(tbl[[paste0(fvars[i], "_2")]], " - ", tbl[[paste0(fvars[i], "_1")]])
      )
    }

    explode <- function(x) {
      m <- do.call(rbind, strsplit(as.character(x), ":", fixed = TRUE))
      out <- as.data.frame(m, stringsAsFactors = FALSE)
      names(out) <- fvars
      out
    }

    g1 <- explode(tbl$Group1)
    g2 <- explode(tbl$Group2)

    same <- g1 == g2
    tbl$facet <- apply(same, 1, function(r) {
      if (all(r))           "all"
      else if (!any(r))     "Effect of other"
      else                  paste0("Effect of ", paste(fvars[!r], collapse = "+"))
    })

    tbl$facet <- .relevel_if_present(tbl$facet, "Effect of other")
    tbl$facet <- factor(tbl$facet, levels = rev(levels(tbl$facet)))

    # simplified labels (your revised version is fine)
    tbl$x_simplified <- sapply(
      seq_len(nrow(tbl)),
      function(i) {
        valid_vars <- fvars[fvars %in% names(g1) & fvars %in% names(g2)]
        if (length(valid_vars) == 0) return("no-vars")

        dif <- character(0); sam <- character(0)
        for (var in valid_vars) {
          g1_val <- g1[i, var]; g2_val <- g2[i, var]
          if (is.na(g1_val) || is.na(g2_val) || length(g1_val) == 0 || length(g2_val) == 0) next
          if (g1_val == g2_val) sam <- c(sam, var) else dif <- c(dif, var)
        }

        if (length(dif) == 0) {
          "no-diff"
        } else if (length(dif) == 1) {
          g1_val <- g1[i, dif]; g2_val <- g2[i, dif]
          if (length(sam) == 0) {
            sprintf("%s: %s vs. %s", dif, g2_val, g1_val)
          } else {
            sam_vals <- sapply(sam, function(s) g2[i, s])
            sam_collapsed <- paste(sprintf("%s: %s", sam, sam_vals), collapse = ", ")
            sprintf("%s\n%s: %s vs. %s", sam_collapsed, dif, g2_val, g1_val)
          }
        } else {
          g1_vals <- sapply(valid_vars, function(var) g1[i, var])
          g2_vals <- sapply(valid_vars, function(var) g2[i, var])
          paste(paste(g1_vals, collapse = ":"), paste(g2_vals, collapse = ":"), sep = " - ")
        }
      },
      USE.NAMES = FALSE
    )

    pd <- ggplot2::position_dodge(width = 0.6)

    p <- ggplot2::ggplot(
      tbl,
      ggplot2::aes(
        x = x_simplified,
        y = log2Ratio_MIC,
        ymin = CI_Lower, ymax = CI_Upper
      )
    )
    if (!is.null(color_by) && color_by %in% names(tbl)) {
      p <- p + ggplot2::aes(color = .data[[color_by]], fill = .data[[color_by]])
    }

    p <- p +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
      ggplot2::geom_col(position = pd, width = 0.5, alpha = 0.35, color = NA) +
      ggplot2::geom_errorbar(position = pd, width = 0.25) +
      ggplot2::geom_point(position = pd, size = 2) +
      ggplot2::coord_flip() +
      ggplot2::facet_grid(rows = ggplot2::vars(facet), scales = "free_y", space = "free") +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        x = "Contrast", y = "MIC ratio (log2 scale)",
        color = color_by, fill = color_by,
        title = "Pairwise log2 MIC* ratios"
      ) +
      ggplot2::theme(
        strip.text.y = ggplot2::element_text(face = "bold"),
        legend.position = if (is.null(color_by)) "none" else "bottom"
      )

    return(p)


  } else if (type == "DoD_ratio") {
    ## ------------------Function to get pairs--------------------------------
    tbl  <- object$dod_ratio_results          # usually one row
    fvars <- names(object$newdata)             # e.g. c("strain","treatment")

    if (length(fvars) < 2L) {
      stop("DoD plots require at least 2 grouping factors.", call. = FALSE)
    }


    if (!"label" %in% names(tbl) &&
        all(c("var1","var2","var1_lvlA","var1_lvlB","var2_lvlC","var2_lvlD") %in% names(tbl))) {
      tbl$label <- sprintf("%s: %s vs %s x %s: %s vs %s",
                           tbl$var1, tbl$var1_lvlB, tbl$var1_lvlA,
                           tbl$var2, tbl$var2_lvlD, tbl$var2_lvlC)
    }

    tbl <- subset(tbl, is.finite(DDlog2MIC) & is.finite(CI_Lower) & is.finite(CI_Upper) #&
                    # DDlog2MIC > 0 & CI_Lower > 0 & CI_Upper > 0
                  )
    if (nrow(tbl) == 0L) stop("No finite DoD ratios to plot.")

    if (nrow(tbl) == 0L) {
      stop("No pairwise difference of difference ratios found. Check that the model has at least two factors.")
    }

    ## ---------------------------------------------------------------
    ## 3.  aesthetic options
    aes_extra <- list(
      fill = tbl$label,
      color = tbl$label
      #fill = "#3182BD", color = "#3182BD"
    )

    ## ---------------------------------------------------------------
    ## 4.  plot
    ggplot2::ggplot(
      tbl,
      ggplot2::aes(x = .data[["label"]],
                   y = DDlog2MIC,
                   ymin = CI_Lower,
                   ymax = CI_Upper,
                   !!!aes_extra)
    ) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                          color = "grey60") +
      ggplot2::geom_errorbar(width = .25, linewidth = .6) +
      ggplot2::geom_point(size = 3) +
      ggplot2::geom_col(width = 0.5, alpha = 0.35, color = NA) +
      ggplot2::coord_flip() +
      ggplot2::theme_minimal() +
      ggplot2::labs(x = NULL,
                    y = "Difference-of-differences (log2 ratio)",
                    title = "DoD on ratio scale") +
      ggplot2::theme(
        axis.text.y  = ggplot2::element_text(face = "bold"),
        legend.position = "none",
        plot.margin  = ggplot2::margin(r = 20)
      )


  } else if (type == "DoD_delta") {

    tbl  <- object$dod_delta_results          # usually one row
    fvars <- names(object$newdata)             # e.g. c("strain","treatment")

    if (length(fvars) < 2L) {
      stop("DoD plots require at least 2 grouping factors.", call. = FALSE)
    }

    if (!"label" %in% names(tbl) &&
        all(c("var1","var2","var1_lvlA","var1_lvlB","var2_lvlC","var2_lvlD") %in% names(tbl))) {
      tbl$label <- sprintf("%s: %s vs %s x %s: %s vs %s",
                           tbl$var1, tbl$var1_lvlB, tbl$var1_lvlA,
                           tbl$var2, tbl$var2_lvlD, tbl$var2_lvlC)
    }

    if (is.null(tbl)) {
      stop("No pairwise difference of differences found. Check that the model has at least two factors.")
    }

    aes_extra <- list(
      fill = tbl$label,
      color = tbl$label
      #fill = "#3182BD", color = "#3182BD"
      )

    ####  plot
    ggplot2::ggplot(
      tbl,
      ggplot2::aes(x = .data[["label"]],
                   y = DDMIC,
                   ymin = CI_Lower,
                   ymax = CI_Upper,
                   !!!aes_extra)
    ) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                          color = "grey60") +
      ggplot2::geom_errorbar(width = .25, linewidth = .6) +
      ggplot2::geom_point(size = 3) +
      ggplot2::geom_col(width = 0.5, alpha = 0.35, color = NA) +
      ggplot2::coord_flip() +
      ggplot2::theme_minimal() +
      ggplot2::labs(x = NULL,
                    y = "Difference-of-differences",
                    title = "DoD of raw MIC") +
      ggplot2::theme(
        axis.text.y  = ggplot2::element_text(face = "bold"),
        legend.position = "none",
        plot.margin  = ggplot2::margin(r = 20)
      )
  }
}

