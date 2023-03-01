plot_dendrogram <-
  function (x, l_balances, type = c("rectangle", "triangle"), center = FALSE,
            edge.root = is.leaf(x) || !is.null(attr(x, "edgetext")),
            nodePar = NULL, edgePar = list(),
            leaflab = c("perpendicular", "textlike", "none"), dLeaf = NULL,
            xlab = "", ylab = "", xaxt="n", yaxt="s",
            horiz = FALSE, frame.plot = FALSE, xlim, ylim, ...)
  {
    type <- match.arg(type)
    leaflab <- match.arg(leaflab)
    hgt <- attr(x, "height")
    if (edge.root && is.logical(edge.root))
      edge.root <- 0.0625 * if(is.leaf(x)) 1 else hgt
    mem.x <- .memberDend(x)
    yTop <- hgt + edge.root
    if(center) { x1 <- 0.5 ; x2 <- mem.x + 0.5 }
    else       { x1 <- 1   ; x2 <- mem.x }
    xl. <- c(x1 - 1/2, x2 + 1/2)
    yl. <- c(0, yTop)
    if (horiz) {## swap and reverse direction on `x':
      tmp <- xl.; xl. <- rev(yl.); yl. <- tmp
      tmp <- xaxt; xaxt <- yaxt; yaxt <- tmp
    }
    if(missing(xlim) || is.null(xlim)) xlim <- xl.
    if(missing(ylim) || is.null(ylim)) ylim <- yl.
    dev.hold(); on.exit(dev.flush())
    plot(0, xlim = xlim, ylim = ylim, type = "n", xlab = xlab, ylab = ylab,
         xaxt = xaxt, yaxt = yaxt, frame.plot = frame.plot, ...)
    if(is.null(dLeaf))
      dLeaf <- .75*(if(horiz) strwidth("w") else strheight("x"))

    if (edge.root) {
      ### FIXME: the first edge + edgetext is drawn here, all others in plotNode()
      ### -----  maybe use trick with adding a single parent node to the top ?
      x0 <- plotNodeLimit(x1, x2, x, center)$x
      if (horiz)
        segments(hgt, x0, yTop, x0)
      else segments(x0, hgt, x0, yTop)
      if (!is.null(et <- attr(x, "edgetext"))) {
        my <- mean(hgt, yTop)
        if (horiz)
          text(my, x0, et)
        else text(x0, my, et)
      }
    }
    plotNode(x1, x2, x, l_balances = l_balances, type = type, center = center, leaflab = leaflab,
             dLeaf = dLeaf, nodePar = nodePar, edgePar = edgePar, horiz = horiz)
  }

### the work horse: plot node (if pch) and lines to all children
plotNode <- function(x1, x2, subtree, l_balances, type, center, leaflab, dLeaf,
                     nodePar, edgePar, horiz = FALSE)
{
  wholetree <- subtree
  depth <- 0L
  llimit <- list()
  KK <- integer()
  kk <- integer()

  repeat {
    inner <- !is.leaf(subtree) && x1 != x2
    yTop <- attr(subtree, "height")
    bx <- plotNodeLimit(x1, x2, subtree, center)
    xTop <- bx$x
    depth <- depth + 1L
    llimit[[depth]] <- bx$limit

    ## handle node specific parameters in "nodePar":
    hasP <- !is.null(nPar <- attr(subtree, "nodePar"))
    if(!hasP) nPar <- nodePar

    if(getOption("verbose")) {
      cat(if(inner)"inner node" else "leaf", ":")
      if(!is.null(nPar)) { cat(" with node pars\n"); str(nPar) }
      cat(if(inner )paste(" height", formatC(yTop),"; "),
          "(x1,x2)= (", formatC(x1, width = 4), ",", formatC(x2, width = 4), ")",
          "--> xTop=", formatC(xTop, width = 8), "\n", sep = "")
    }

    Xtract <- function(nam, L, default, indx)
      rep(if(nam %in% names(L)) L[[nam]] else default,
          length.out = indx)[indx]
    asTxt <- function(x) # to allow 'plotmath' labels:
      if(is.character(x) || is.expression(x) || is.null(x)) x else as.character(x)

    i <- if(inner || hasP) 1 else 2 # only 1 node specific par

    if(!is.null(nPar)) { ## draw this node
      pch <- Xtract("pch", nPar, default = 1L:2,	 i)
      cex <- Xtract("cex", nPar, default = c(1,1),	 i)
      col <- Xtract("col", nPar, default = par("col"), i)
      bg <- Xtract("bg", nPar, default = par("bg"), i)
      points(if (horiz) cbind(yTop, xTop) else cbind(xTop, yTop),
             pch = pch, bg = bg, col = col, cex = cex)
    }

    if(leaflab == "textlike")
      p.col <- Xtract("p.col", nPar, default = "white", i)
    lab.col <- Xtract("lab.col", nPar, default = par("col"), i)
    lab.cex <- Xtract("lab.cex", nPar, default = c(1,1), i)
    lab.font <- Xtract("lab.font", nPar, default = par("font"), i)
    lab.xpd <- Xtract("xpd", nPar, default = c(TRUE, TRUE), i)
    if (is.leaf(subtree)) {
      ## label leaf
      if (leaflab == "perpendicular") { # somewhat like plot.hclust
        if(horiz) {
          X <- yTop + dLeaf * lab.cex
          Y <- xTop; srt <- 0; adj <- c(0, 0.5)
        }
        else {
          Y <- yTop - dLeaf * lab.cex
          X <- xTop; srt <- 90; adj <- 1
        }
        nodeText <- asTxt(attr(subtree,"label"))
        text(X, Y, nodeText, xpd = lab.xpd, srt = srt, adj = adj,
             cex = lab.cex, col = lab.col, font = lab.font)
      }
    }
    else if (inner) {
      segmentsHV <- function(x0, y0, x1, y1) {
        if (horiz)
          segments(y0, x0, y1, x1, col = col, lty = lty, lwd = lwd)
        else segments(x0, y0, x1, y1, col = col, lty = lty, lwd = lwd)
      }
      for (k in seq_along(subtree)) {
        child <- subtree[[k]]
        ## draw lines to the children and draw them recursively
        yBot <- attr(child, "height")
        if (getOption("verbose")) cat("ch.", k, "@ h=", yBot, "; ")
        if (is.null(yBot))
          yBot <- 0
        xBot <-
          if (center) mean(bx$limit[k:(k + 1)])
        else bx$limit[k] + .midDend(child)

        hasE <- !is.null(ePar <- attr(child, "edgePar"))
        if (!hasE)
          ePar <- edgePar
        i <- if (!is.leaf(child) || hasE) 1 else 2
        ## define line attributes for segmentsHV():
        col <- Xtract("col", ePar, default = par("col"), i)
        lty <- Xtract("lty", ePar, default = par("lty"), i)
        lwd <- Xtract("lwd", ePar, default = par("lwd"), i)
        if (type == "triangle") {
          segmentsHV(xTop, yTop, xBot, yBot)
        }
        else { # rectangle
          segmentsHV(xTop,yTop, xBot,yTop)# h
          segmentsHV(xBot,yTop, xBot,yBot)# v
        }
        vln <- NULL
        if (is.leaf(child) && leaflab == "textlike") {
          nodeText <- asTxt(attr(child,"label"))
          if(getOption("verbose"))
            cat('-- with "label"',format(nodeText))
          hln <- 0.6 * strwidth(nodeText, cex = lab.cex)/2
          vln <- 1.5 * strheight(nodeText, cex = lab.cex)/2
          rect(xBot - hln, yBot,
               xBot + hln, yBot + 2 * vln, col = p.col)
          text(xBot, yBot + vln, nodeText, xpd = lab.xpd,
               cex = lab.cex, col = lab.col, font = lab.font)
        }
        if (!is.null(attr(child, "edgetext"))) {
          edgeText <- asTxt(attr(child, "edgetext"))
          if(getOption("verbose"))
            cat('-- with "edgetext"',format(edgeText))
          if (!is.null(vln)) {
            mx <-
              if(type == "triangle")
                (xTop+ xBot+ ((xTop - xBot)/(yTop - yBot)) * vln)/2
            else xBot
            my <- (yTop + yBot + 2 * vln)/2
          }
          else {
            mx <- if(type == "triangle") (xTop + xBot)/2 else xBot
            my <- (yTop + yBot)/2
          }
          ## Both for "triangle" and "rectangle" : Diamond + Text

          p.col <- Xtract("p.col", ePar, default = "white", i)
          p.border <- Xtract("p.border", ePar, default = par("fg"), i)
          ## edge label pars: defaults from the segments pars
          p.lwd <- Xtract("p.lwd", ePar, default = lwd, i)
          p.lty <- Xtract("p.lty", ePar, default = lty, i)
          t.col <- Xtract("t.col", ePar, default = col, i)
          t.cex <- Xtract("t.cex", ePar, default =  1,  i)
          t.font <- Xtract("t.font", ePar, default = par("font"), i)

          vlm <- strheight(c(edgeText,"h"), cex = t.cex)/2
          hlm <- strwidth (c(edgeText,"m"), cex = t.cex)/2
          hl3 <- c(hlm[1L], hlm[1L] + hlm[2L], hlm[1L])
          if(horiz) {
            polygon(my+ c(-hl3, hl3), mx + sum(vlm)*c(-1L:1L, 1L:-1L),
                    col = p.col, border = p.border,
                    lty = p.lty, lwd = p.lwd)
            text(my, mx, edgeText, cex = t.cex, col = t.col,
                 font = t.font)
          } else {
            polygon(mx+ c(-hl3, hl3), my + sum(vlm)*c(-1L:1L, 1L:-1L),
                    col = p.col, border = p.border,
                    lty = p.lty, lwd = p.lwd)
            text(mx, my, edgeText, cex = t.cex, col = t.col,
                 font = t.font)
          }
        }
      }
      labels_ = unlist(dendrapply(subtree, function(x) if(is.leaf(x)) return(attr(x, "label")) ))
      ibal = which(sapply(l_balances, function(x) length(x) == length(labels_) & all(x %in% labels_)))
      points(xTop, yTop, col = 'white', cex = 5, pch = 19)
      text(xTop, yTop, names(l_balances)[ibal])
    }

    if (inner && length(subtree)) {
      KK[depth] <- length(subtree)
      if (storage.mode(kk) != storage.mode(KK))
        storage.mode(kk) <- storage.mode(KK)

      ## go to first child
      kk[depth] <- 1L
      x1 <- bx$limit[1L]
      x2 <- bx$limit[2L]
      subtree <- subtree[[1L]]
    }
    else {
      repeat {
        depth <- depth - 1L
        if (!depth || kk[depth] < KK[depth]) break
      }
      if (!depth) break
      length(kk) <- depth
      kk[depth] <- k <- kk[depth] + 1L
      x1 <- llimit[[depth]][k]
      x2 <- llimit[[depth]][k + 1L]
      subtree <- wholetree[[kk]]
    }
  } ## repeat
  invisible()
}

plotNodeLimit <- function(x1, x2, subtree, center)
{
  ## get the left borders limit[k] of all children k=1..K, and
  ## the handle point `x' for the edge connecting to the parent.
  inner <- !is.leaf(subtree) && x1 != x2
  limit <- c(x1,
             if(inner) {
               K <- length(subtree)
               mTop <- .memberDend(subtree)
               limit <- integer(K)
               xx1 <- x1
               for(k in 1L:K) {
                 m <- .memberDend(subtree[[k]])
                 ##if(is.null(m)) m <- 1
                 xx1 <- xx1 + (if(center) (x2-x1) * m/mTop else m)
                 limit[k] <- xx1
               }
               limit
             } else ## leaf
               x2)
  mid <- attr(subtree, "midpoint")
  center <- center || (inner && !is.numeric(mid))
  x <- if(center) mean(c(x1,x2)) else x1 + (if(inner) mid else 0)
  list(x = x, limit = limit)
}
.memberDend <- function(x) {
  r <- attr(x,"x.member")
  if(is.null(r)) {
    r <- attr(x,"members")
    if(is.null(r)) r <- 1L
  }
  r
}
.midDend <- function(x)
  if(is.null(mp <- attr(x, "midpoint"))) 0 else mp
