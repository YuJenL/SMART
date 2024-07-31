osIsWin <- Sys.info()["sysname"]=="Windows"
isIcon <- is.element("icon",ls())

peakabunCheck <- function(rawdata){
 info <- rawdata[,1:3]
 newdata <- as.matrix(rawdata[,-(1:3)])

 
 idx.i <- which(!colSds(newdata,na.rm=T)%in%c(0,NA))
 newdata <- newdata[,idx.i]
 
 idx.p <- which(!rowSds(newdata,na.rm=T)%in%c(0,NA))
 newdata <- cbind(info[idx.p,],newdata[idx.p,])
 
 newdata
}

Meta_3D.plot = function(meta_mass = mass, meta_peak = peak, meta_RT = RT, rt_begin = 20, rt_range = 50, m_reso=100, t_reso = 100, u_inten = Inf, l_inten = -Inf, h_angle = 15, v_angle = 15, mz_begin = 100, mz_range = 100){
	
	scan_idx = which(meta_RT>rt_begin & meta_RT<(rt_begin+rt_range))
	meta_mass = meta_mass[scan_idx]
	meta_peak = meta_peak[scan_idx]
	meta_mass_idx = sapply(meta_mass, function(x) which(x>mz_begin & x<(mz_begin+mz_range)))
	meta_mass = sapply(1:length(scan_idx), function(x) meta_mass[[x]][meta_mass_idx[[x]]])
	meta_peak = sapply(1:length(scan_idx), function(x) meta_peak[[x]][meta_mass_idx[[x]]])
	meta_RT = meta_RT[scan_idx]
	mass_range = sapply(meta_mass, range)
	mass_fix = seq(min(mass_range), max(mass_range), len = m_reso)
	mass_size = mass_fix[2]-mass_fix[1]
	mass_begin = mass_fix[1]-mass_size/2
	mass_idx = sapply(meta_mass, function(x) ceiling((x-mass_begin)/mass_size))
	peak_max = matrix(0, m_reso, length(meta_RT))
	for(i in 1:length(meta_RT)){
		temp = as.integer(names(table(mass_idx[[i]])))
		peak_max[temp,i] = sapply(temp, function(x) max(meta_peak[[i]][which(mass_idx[[i]]==x)]))
	}
	RT_fix = seq(min(meta_RT), max(meta_RT), len = t_reso)
	RT_size = RT_fix[2]-RT_fix[1]
	RT_begin = RT_fix[1]-RT_size/2
	RT_idx = ceiling((meta_RT-RT_begin)/RT_size)
	peak_max_RT = matrix(0, m_reso, t_reso)
	for(i in 1:m_reso){
		peak_max_RT[i,] = sapply(1:t_reso, function(x) ifelse(length(which(RT_idx==x))==0, 0, max(peak_max[i,which(RT_idx==x)])))
	}
	peak_max_RT[peak_max_RT>u_inten] = min(peak_max_RT)
	peak_max_RT[peak_max_RT<l_inten] = l_inten
	if(is.null(dev.list)) switch(ws, par(mar = c(3, 2, -0.1, 1)+0.1, pin = c(570/96,620/96)), par(mar = c(5, 4, 2, 2)+0.1, pin = c(540/96,410/96)), par(mar = c(5, 4, 2, 2)+0.1, pin = c(510/96,367/96))) 
	persp(RT_fix, mass_fix, t(peak_max_RT), border = NA, col = "lawngreen", theta = h_angle, phi = v_angle, shade = 0.75, ticktype= "detailed", xlab = "Retention Time", ylab = "m/z", zlab = "Intensity", main = NAME, expand = 0.7)
	
}

show3D <- function(){
	
	tkconfigure(tt,cursor="watch")
    a = openMSfile(tclvalue(tcl(treeWidget,"selection","get")))
	NAME <<- gsub(".*/(.*?).mzXML","\\1", tclvalue(tcl(treeWidget,"selection","get")))
	NAME <<- gsub("\\d{4,8}_(*?)", "\\1", NAME)
	hd = header(a)
	RT <<- hd[,"retentionTime"]
	a_peak = peaks(a)
	peak <<- sapply(a_peak, function(x) x[,2])
	mass <<- sapply(a_peak, function(x) x[,1])

	ocl <<- cl <<- substitute(Meta_3D.plot)
    exargs <<- as.list(quote(list()))
    replot <- function() {eval(cl)}

	if(exists("img")) tkpack.forget(img)
    img <<- tkrplot(frameforplot, replot, hscale = 2, vscale = 1.75)
    tkpack(img, side = "top")
	
	tkconfigure(tl_m, state = "normal")
	tkconfigure(tl_t, state = "normal")
	tkconfigure(tl_ub, state = "normal")
	tkconfigure(tl_lb, state = "normal")
	
	if(tclvalue(tkwinfo("exists", ts_RT_begin))=="1") tkdestroy(ts_RT_begin)
	if(tclvalue(tkwinfo("exists", ts_RT_begin_23))=="1") tkpack(ts_RT_begin_23)

	if(tclvalue(tkwinfo("exists", ts_RT_begin_23))=="0"){
		tkpack(ts_RT_begin_23 <<- tkscale(fr3.1, width = "14", variable = RT_begin, orient = "horizontal"))
	}
	tkconfigure(ts_RT_begin_23, from = 0, to = min(max(RT))-0.5, resolution = 0.1, state = "normal")
	tkbind(ts_RT_begin_23, "<ButtonRelease-1>", function(...) tkrreplot(img, hscale = switch(ws, as.numeric(tclvalue(hsc))-0.2, as.numeric(tclvalue(hsc))-0.35, as.numeric(tclvalue(hsc))-0.5), vscale = switch(ws, as.numeric(tclvalue(vsc)), as.numeric(tclvalue(vsc))-0.5, as.numeric(tclvalue(vsc))-0.5)))
	
	
	if(tclvalue(tkwinfo("exists", ts_RT_range))=="1") tkdestroy(ts_RT_range)
	if(tclvalue(tkwinfo("exists", ts_RT_range_23))=="1") tkpack(ts_RT_range_23)

	if(tclvalue(tkwinfo("exists", ts_RT_range_23))=="0"){
		tkpack(ts_RT_range_23 <<- tkscale(fr3.2, width = "14", variable = RT_range, orient = "horizontal"))
	}
	tkconfigure(ts_RT_range_23, from = 2, to = min(max(RT)), resolution = 0.1, state = "normal")
	tkbind(ts_RT_range_23, "<ButtonRelease-1>", function(...) tkrreplot(img, hscale = switch(ws, as.numeric(tclvalue(hsc))-0.2, as.numeric(tclvalue(hsc))-0.35, as.numeric(tclvalue(hsc))-0.5), vscale = switch(ws, as.numeric(tclvalue(vsc)), as.numeric(tclvalue(vsc))-0.5, as.numeric(tclvalue(vsc))-0.5)))
	
	
	tkconfigure(ts_MZ_begin, from = 0, to = max(sapply(mass, max))-1, resolution = 0.1, state = "normal")
	tkbind(ts_MZ_begin, "<ButtonRelease-1>", function(...) tkrreplot(img, hscale = switch(ws, as.numeric(tclvalue(hsc))-0.2, as.numeric(tclvalue(hsc))-0.35, as.numeric(tclvalue(hsc))-0.5), vscale = switch(ws, as.numeric(tclvalue(vsc)), as.numeric(tclvalue(vsc))-0.5, as.numeric(tclvalue(vsc))-0.5)))
	
	
	tkconfigure(ts_MZ_range, from = 2, to = max(sapply(mass, max)), resolution = 0.1, state = "normal")
	tkbind(ts_MZ_range, "<ButtonRelease-1>", function(...) tkrreplot(img, hscale = switch(ws, as.numeric(tclvalue(hsc))-0.2, as.numeric(tclvalue(hsc))-0.35, as.numeric(tclvalue(hsc))-0.5), vscale = switch(ws, as.numeric(tclvalue(vsc)), as.numeric(tclvalue(vsc))-0.5, as.numeric(tclvalue(vsc))-0.5)))
	
	
	tkconfigure(ts_v_angle, from = 0, to = 90, resolution = 1, state = "normal")
	tkbind(ts_v_angle, "<ButtonRelease-1>", function(...) tkrreplot(img, hscale = switch(ws, as.numeric(tclvalue(hsc))-0.2, as.numeric(tclvalue(hsc))-0.35, as.numeric(tclvalue(hsc))-0.5), vscale = switch(ws, as.numeric(tclvalue(vsc)), as.numeric(tclvalue(vsc))-0.5, as.numeric(tclvalue(vsc))-0.5)))
	
	
	tkconfigure(ts_h_angle, from = 0, to = 90, resolution = 1, state = "normal")
	tkbind(ts_h_angle, "<ButtonRelease-1>", function(...) tkrreplot(img, hscale = switch(ws, as.numeric(tclvalue(hsc))-0.2, as.numeric(tclvalue(hsc))-0.35, as.numeric(tclvalue(hsc))-0.5), vscale = switch(ws, as.numeric(tclvalue(vsc)), as.numeric(tclvalue(vsc))-0.5, as.numeric(tclvalue(vsc))-0.5)))
	
	tmpcl <- as.list(cl)
	tmpl_rt_begin <- list(substitute(as.numeric(tclvalue(vname)), list(vname = as.character(RT_begin))))
	tmpl_rt_range <- list(substitute(as.numeric(tclvalue(vname)), list(vname = as.character(RT_range))))
	tmpl_mz_begin <- list(substitute(as.numeric(tclvalue(vname)), list(vname = as.character(MZ_begin))))
	tmpl_mz_range <- list(substitute(as.numeric(tclvalue(vname)), list(vname = as.character(MZ_range))))
	tmpl_v_angle <- list(substitute(as.numeric(tclvalue(vname)), list(vname = as.character(v_an))))
	tmpl_h_angle <- list(substitute(as.numeric(tclvalue(vname)), list(vname = as.character(h_an))))
	tmpl_m <- list(substitute(as.numeric(tclvalue(vname)),list(vname = as.character(m_r))))
	tmpl_t <- list(substitute(as.numeric(tclvalue(vname)),list(vname = as.character(t_r))))
	tmpl_lb <- list(substitute(as.numeric(tclvalue(vname)),list(vname = as.character(lowerbound))))
	tmpl_ub <- list(substitute(as.numeric(tclvalue(vname)),list(vname = as.character(upperbound))))
	names(tmpl_m) <- "m_reso"
	names(tmpl_t) <- "t_reso"
	names(tmpl_lb) <- "l_inten"
	names(tmpl_ub) <- "u_inten"
	names(tmpl_rt_begin) <- "rt_begin"
	names(tmpl_rt_range) <- "rt_range"
	names(tmpl_mz_begin) <- "mz_begin"
	names(tmpl_mz_range) <- "mz_range"
	names(tmpl_v_angle) <- "v_angle"
	names(tmpl_h_angle) <- "h_angle"
	cl <<- as.call(c(tmpcl, tmpl_rt_begin, tmpl_rt_range, tmpl_mz_begin, tmpl_mz_range, tmpl_v_angle, tmpl_h_angle, tmpl_m, tmpl_t, tmpl_lb, tmpl_ub))
	exargs <<- c(exargs, tmpl_rt_begin, tmpl_rt_range, tmpl_mz_begin, tmpl_mz_range, tmpl_v_angle, tmpl_h_angle, tmpl_m, tmpl_t, tmpl_lb, tmpl_ub)
	
	
    tkrreplot(img, hscale = switch(ws, as.numeric(tclvalue(hsc))-0.2, as.numeric(tclvalue(hsc))-0.35, as.numeric(tclvalue(hsc))-0.5), vscale = switch(ws, as.numeric(tclvalue(vsc)), as.numeric(tclvalue(vsc))-0.5, as.numeric(tclvalue(vsc))-0.5))
	tkconfigure(refresh.but, command = function() {
		tkrreplot(img, hscale = switch(ws, as.numeric(tclvalue(hsc))-0.2, as.numeric(tclvalue(hsc))-0.35, as.numeric(tclvalue(hsc))-0.5), vscale = switch(ws, as.numeric(tclvalue(vsc)), as.numeric(tclvalue(vsc))-0.5, as.numeric(tclvalue(vsc))-0.5))
	})
	tkconfigure(tt,cursor="arrow")
}

Meta_2D.plot = function(meta_mass = mass, meta_peak = peak, meta_RT = RT, rt_begin = 20, rt_range = 50, m_reso=100, t_reso = 100, u_inten = Inf, l_inten = -Inf, mz_begin = 100, mz_range = 100){
	
	scan_idx = which(meta_RT>rt_begin & meta_RT<(rt_begin+rt_range))
	meta_mass = meta_mass[scan_idx]
	meta_peak = meta_peak[scan_idx]
	meta_mass_idx = sapply(meta_mass, function(x) which(x>mz_begin & x<(mz_begin+mz_range)))
	meta_mass = sapply(1:length(scan_idx), function(x) meta_mass[[x]][meta_mass_idx[[x]]])
	meta_peak = sapply(1:length(scan_idx), function(x) meta_peak[[x]][meta_mass_idx[[x]]])
	meta_RT = meta_RT[scan_idx]
	mass_range = sapply(meta_mass, range)
	mass_fix <<- seq(min(mass_range), max(mass_range), len = m_reso)
	mass_size = mass_fix[2]-mass_fix[1]
	mass_begin = mass_fix[1]-mass_size/2
	mass_idx = sapply(meta_mass, function(x) ceiling((x-mass_begin)/mass_size))
	peak_max = matrix(0, m_reso, length(meta_RT))
	for(i in 1:length(meta_RT)){
		temp = as.integer(names(table(mass_idx[[i]])))
		peak_max[temp,i] = sapply(temp, function(x) max(meta_peak[[i]][which(mass_idx[[i]]==x)]))
	}
	RT_fix <<- seq(min(meta_RT), max(meta_RT), len = t_reso)
	RT_size = RT_fix[2]-RT_fix[1]
	RT_begin = RT_fix[1]-RT_size/2
	RT_idx = ceiling((meta_RT-RT_begin)/RT_size)
	peak_max_RT <<- matrix(0, m_reso, t_reso)
	for(i in 1:m_reso){
		peak_max_RT[i,] <<- sapply(1:t_reso, function(x) ifelse(length(which(RT_idx==x))==0, 0, max(peak_max[i,which(RT_idx==x)])))
	}
	
	peak_max_RT[peak_max_RT>u_inten] = min(peak_max_RT)
	peak_max_RT[peak_max_RT<l_inten] = l_inten
	brk <<- quantile(unique(c(ceiling(peak_max_RT))), probs = seq(0,1,by = 0.05))
	peak_max_RT_ = matrix(as.numeric(cut(peak_max_RT, brk, include.lowest = T, right = F)), nrow = nrow(peak_max_RT))
	
	
	
	
	
	
	if(is.null(dev.list)) switch(ws, par(pin = c(720/96, 600/96)), par(pin = c(540/96, 400/96)), par(pin = c(480/96, 360/96)))
	filled.contour(RT_fix, mass_fix, t(peak_max_RT_), levels = as.numeric(factor(brk)), color.palette = colorRampPalette(c("white", "green3", "red"), space = "rgb"), plot.title = title(main = NAME, xlab = "Retention Time", ylab = "m/z"), key.axes = axis(4, at = seq(1, length(brk), by = 2), labels = round(brk[seq(1, length(brk), by = 2)], 1), cex.lab = 1))


	
}

show2D <- function(){
	tkconfigure(tt,cursor="watch")
	a = openMSfile(tclvalue(tcl(treeWidget,"selection","get")))
	
	NAME <<- gsub(".*/(.*?).mzXML","\\1", tclvalue(tcl(treeWidget,"selection","get")))
	NAME <<- gsub("\\d{4,8}_(*?)", "\\1", NAME)
	hd = header(a)
	RT <<- hd[,"retentionTime"]
	a_peak = peaks(a)
	peak <<- sapply(a_peak, function(x) x[,2])
	mass <<- sapply(a_peak, function(x) x[,1])

	ocl <<- cl <<- substitute(Meta_2D.plot)
    exargs <<- as.list(quote(list()))
    replot <- function() {eval(cl)}
	if(exists("img")) tkpack.forget(img)
    img <<- tkrplot(frameforplot, replot, hscale = 2, vscale = switch(ws, 1.75-0.1, 1.5, 1.4))
    tkpack(img, side = "top")


	tkconfigure(tl_m, state = "normal")
	tkconfigure(tl_t, state = "normal")
	tkconfigure(tl_ub, state = "normal")
	tkconfigure(tl_lb, state = "normal")

	
	if(tclvalue(tkwinfo("exists", ts_RT_begin))=="1") tkdestroy(ts_RT_begin)
	if(tclvalue(tkwinfo("exists", ts_RT_begin_23))=="1") tkpack(ts_RT_begin_23)

	if(tclvalue(tkwinfo("exists", ts_RT_begin_23))=="0"){
		tkpack(ts_RT_begin_23 <<- tkscale(fr3.1, width = "14", variable = RT_begin, orient = "horizontal"))
	}
	tkconfigure(ts_RT_begin_23, from = 0, to = min(max(RT))-0.5, resolution = 0.1, state = "normal")
	tkbind(ts_RT_begin_23, "<ButtonRelease-1>", function(...) tkrreplot(img, hscale = as.numeric(tclvalue(hsc)), vscale = switch(ws, as.numeric(tclvalue(vsc))-0.1, as.numeric(tclvalue(vsc))-0.39, as.numeric(tclvalue(vsc))-0.5)))
	
	
	if(tclvalue(tkwinfo("exists", ts_RT_range))=="1") tkdestroy(ts_RT_range)
	if(tclvalue(tkwinfo("exists", ts_RT_range_23))=="1") tkpack(ts_RT_range_23)

	if(tclvalue(tkwinfo("exists", ts_RT_range_23))=="0"){
		tkpack(ts_RT_range_23 <<- tkscale(fr3.2, width = "14", variable = RT_range, orient = "horizontal"))
	}
	tkconfigure(ts_RT_range_23, from = 2, to = min(max(RT)), resolution = 0.1, state = "normal")
	tkbind(ts_RT_range_23, "<ButtonRelease-1>", function(...) tkrreplot(img, hscale = as.numeric(tclvalue(hsc)), vscale = switch(ws, as.numeric(tclvalue(vsc))-0.1, as.numeric(tclvalue(vsc))-0.39, as.numeric(tclvalue(vsc))-0.5)))
	
	
	tkconfigure(ts_MZ_begin, from = 0, to = max(sapply(mass, max))-1, resolution = 0.1, state = "normal")
	tkbind(ts_MZ_begin, "<ButtonRelease-1>", function(...) tkrreplot(img, hscale = as.numeric(tclvalue(hsc)), vscale = switch(ws, as.numeric(tclvalue(vsc))-0.1, as.numeric(tclvalue(vsc))-0.39, as.numeric(tclvalue(vsc))-0.5)))
	
	
	tkconfigure(ts_MZ_range, from = 2, to = max(sapply(mass, max)), resolution = 0.1, state = "normal")
	tkbind(ts_MZ_range, "<ButtonRelease-1>", function(...) tkrreplot(img, hscale = as.numeric(tclvalue(hsc)), vscale = switch(ws, as.numeric(tclvalue(vsc))-0.1, as.numeric(tclvalue(vsc))-0.39, as.numeric(tclvalue(vsc))-0.5)))
	
	tkconfigure(ts_v_angle, state = "disable")
	tkconfigure(ts_h_angle, state = "disable")
	
	tmpcl <- as.list(cl)
	tmpl_rt_begin <- list(substitute(as.numeric(tclvalue(vname)), list(vname = as.character(RT_begin))))
	tmpl_rt_range <- list(substitute(as.numeric(tclvalue(vname)), list(vname = as.character(RT_range))))
	tmpl_mz_begin <- list(substitute(as.numeric(tclvalue(vname)), list(vname = as.character(MZ_begin))))
	tmpl_mz_range <- list(substitute(as.numeric(tclvalue(vname)), list(vname = as.character(MZ_range))))
	tmpl_m <- list(substitute(as.numeric(tclvalue(vname)),list(vname = as.character(m_r))))
	tmpl_t <- list(substitute(as.numeric(tclvalue(vname)),list(vname = as.character(t_r))))
	tmpl_lb <- list(substitute(as.numeric(tclvalue(vname)),list(vname = as.character(lowerbound))))
	tmpl_ub <- list(substitute(as.numeric(tclvalue(vname)),list(vname = as.character(upperbound))))

	names(tmpl_m) <- "m_reso"
	names(tmpl_t) <- "t_reso"
	names(tmpl_lb) <- "l_inten"
	names(tmpl_ub) <- "u_inten"
	names(tmpl_rt_begin) <- "rt_begin"
	names(tmpl_rt_range) <- "rt_range"
	names(tmpl_mz_begin) <- "mz_begin"
	names(tmpl_mz_range) <- "mz_range"
	cl <<- as.call(c(tmpcl, tmpl_rt_begin, tmpl_rt_range, tmpl_mz_begin, tmpl_mz_range, tmpl_m, tmpl_t, tmpl_lb, tmpl_ub))
	exargs <<- c(exargs, tmpl_rt_begin, tmpl_rt_range, tmpl_mz_begin, tmpl_mz_range, tmpl_m, tmpl_t, tmpl_lb, tmpl_ub)

    tkrreplot(img, hscale = as.numeric(tclvalue(hsc)), vscale = switch(ws, as.numeric(tclvalue(vsc))-0.1, as.numeric(tclvalue(vsc))-0.39, as.numeric(tclvalue(vsc))-0.5))
	tkconfigure(refresh.but, command = function() {
		tkrreplot(img, hscale = as.numeric(tclvalue(hsc)), vscale = switch(ws, as.numeric(tclvalue(vsc))-0.1, as.numeric(tclvalue(vsc))-0.39, as.numeric(tclvalue(vsc))-0.5))
	})
	tkconfigure(tt,cursor="arrow")
}


TICplot = function(meta_NAME = NAME, meta_mass = mass, meta_peak = peak, meta_RT = RT, meta_TIC = TIC, rt_begin = 20, rt_range = 50, cex = 0.8, inset = c(-0.2, 0.08), cex.axis = 1, cex.main = 1.2, cex.lab = 1){
	totalRT <<- NULL
	totalTIC <<- NULL
	if(length(meta_NAME)==1){
		RT_idx <<- which(meta_RT>rt_begin & meta_RT<(rt_begin+rt_range))
		totalRT <<- meta_RT[RT_idx]
		ind_mass <<- meta_mass[RT_idx]
		ind_peak <<- meta_peak[RT_idx]
		totalTIC <<- as.data.frame(meta_TIC[RT_idx])
	}else{
		RT_idx <<- which(meta_RT>rt_begin & meta_RT<(rt_begin+rt_range))
		totalRT <<- meta_RT[RT_idx]
		totalTIC <<- meta_TIC[RT_idx,]
	}
	
		
		
		
		
		
		
		
		
		
		
		
		
	
	
	if(is.null(dev.list)){
		switch(ws, par(mar = c(5+1, 4+5.5, 2, 2), pin = c(540/96,500/96)), par(mar = c(5-1, 4+3, 2, 2), pin = c(540/96,367/96)), par(mar = c(5-1, 4+3, 2, 2), pin = c(460/96,320/96)))
	}else{
		par(mar = c(5-1, 4+5.5, 4+2, 2))
	}

	plot(x = c(rt_begin, rt_begin + rt_range), y = c(0, max(sapply(totalTIC, max))), type = "n", xlim = c(rt_begin, rt_begin + rt_range), ylim = c(0, max(sapply(totalTIC, max))), xlab = "Retention Time (sec)", ylab = "Total Ion Current", main = "TIC plot", cex.axis = cex.axis, cex.main = cex.main, cex.lab = cex.lab)
	par(xpd = NA)
	
	legend(x = "bottomleft", legend = meta_NAME, col = rainbow(length(meta_NAME)), pch = 16, bty = "o", box.col = "gray", cex = cex, inset = inset)
	par(new = T)
	for(f in 1:length(meta_NAME)) lines(totalRT, totalTIC[,f], col = rainbow(length(meta_NAME))[f])
	parPlotSize <<- par("plt")
	usrCoords	<<- par("usr")
}

labelClosestPoint <- function(xClick,yClick,imgXcoords,imgYcoords)
{
	squared.Distance <- (xClick-imgXcoords)^2 + (yClick-imgYcoords)^2
	indexClosest <- which.min(squared.Distance)
	if(exists("plotpeak")) tkdestroy(plotpeak)
	plotpeak <<- tktoplevel(); if(isIcon) tk2ico.set(plotpeak,icon)
	tkwm.title(plotpeak, "Spectra plot")
	replot <- function() plot(mass[[RT_idx[indexClosest]]], peak[[RT_idx[indexClosest]]], type = "l", xlab = "m/z", ylab = "Intensity", main = paste(NAME, "   Scan #", RT_idx[indexClosest], "    RT ", totalRT[indexClosest], " secs", sep = ""))
	imgpeak <<- tkrplot(plotpeak, replot, hscale = 2, vscale = 1.75)
	tkpack(imgpeak, side = "top")
}

OnLeftClick <- function(x,y)
{
	xClick <- x
	yClick <- y

	width <- as.numeric(tclvalue(tkwinfo("reqwidth",img)))
	height <- as.numeric(tclvalue(tkwinfo("reqheight",img)))

	xMin <- parPlotSize[1] * width
	xMax <- parPlotSize[2] * width
	yMin <- parPlotSize[3] * height
	yMax <- parPlotSize[4] * height

	rangeX <- usrCoords[2] - usrCoords[1]
	rangeY <- usrCoords[4] - usrCoords[3]

	imgXcoords <- (totalRT-usrCoords[1])*(xMax-xMin)/rangeX + xMin
	imgYcoords <- (totalTIC-usrCoords[3])*(yMax-yMin)/rangeY + yMin

	xClick <- as.numeric(xClick)+0.5
	yClick <- as.numeric(yClick)+0.5
	yClick <- height - yClick

	xPlotCoord <- usrCoords[1]+(xClick-xMin)*rangeX/(xMax-xMin)
	yPlotCoord <- usrCoords[3]+(yClick-yMin)*rangeY/(yMax-yMin)

	labelClosestPoint(xClick,yClick,imgXcoords,imgYcoords)
}

showTIC <- function(){
	tkconfigure(tt,cursor="watch")
	FILE <- unlist(strsplit(tclvalue(tcl(treeWidget,"selection","get")), split = " |,,"))
	NAME <<- gsub(".*/(.*?).mzXML","\\1", FILE)
	NAME <<- gsub("\\d{4,8}_(*?)", "\\1", NAME)
	if(length(NAME)==1){
		a = openMSfile(FILE)
		hd = header(a)
		RT <<- hd[,"retentionTime"]
		TIC <<- hd[,"totIonCurrent"]
		a_peak = peaks(a)
		peak <<- sapply(a_peak, function(x) x[,2])
		mass <<- sapply(a_peak, function(x) x[,1])
	}else{
		TIC <<- NULL
		RT <<- NULL
		for(i in 1:length(NAME)){
			a = openMSfile(FILE[i])
			hd = header(a)
			RT <<- cbind(RT, hd[,"retentionTime"])
			TIC <<- cbind(TIC, hd[,"totIonCurrent"])
		}
		RT <<- rowMeans(RT)
	}
	
	ocl <<- cl <<- substitute(TICplot)
    exargs <<- as.list(quote(list()))
    replot <- function() {eval(cl)}
	if(exists("img")) tkpack.forget(img)
    img <<- tkrplot(frameforplot, replot, hscale = 2, vscale = 1.75)
    tkpack(img, side = "top")


	tkconfigure(tl_m, state = "disable")
	tkconfigure(tl_t, state = "disable")
	tkconfigure(tl_lb, state = "disable")
	tkconfigure(tl_ub, state = "disable")

	
	
	if(tclvalue(tkwinfo("exists", ts_RT_begin_23))=="1") tkdestroy(ts_RT_begin_23)
	if(tclvalue(tkwinfo("exists", ts_RT_begin))=="1") tkpack(ts_RT_begin)

	if(tclvalue(tkwinfo("exists", ts_RT_begin))=="0"){
		tkpack(ts_RT_begin <<- tkscale(fr3.1, width = "14", variable = RT_begin, orient = "horizontal"))
	}
	tkconfigure(ts_RT_begin, from = 0, to = min(max(RT))-0.5, resolution = 0.1, state = "normal")
	
	tkbind(ts_RT_begin, "<B1-Motion>", function(...) tkrreplot(img, hscale = switch(ws, as.numeric(tclvalue(hsc)), as.numeric(tclvalue(hsc)), as.numeric(tclvalue(hsc))-0.1), vscale = switch(ws, as.numeric(tclvalue(vsc))-0.08, as.numeric(tclvalue(vsc))-0.4, as.numeric(tclvalue(vsc))-0.4)))
	
	
	if(tclvalue(tkwinfo("exists", ts_RT_range_23))=="1") tkdestroy(ts_RT_range_23)
	if(tclvalue(tkwinfo("exists", ts_RT_range))=="1") tkpack(ts_RT_range)

	if(tclvalue(tkwinfo("exists", ts_RT_range))=="0"){
		tkpack(ts_RT_range <<- tkscale(fr3.2, width = "14", variable = RT_range, orient = "horizontal"))
	}
	tkconfigure(ts_RT_range, from = 2, to = min(max(RT)), resolution = 0.1, state = "normal")
	tkbind(ts_RT_range, "<B1-Motion>", function(...) tkrreplot(img, hscale = switch(ws, as.numeric(tclvalue(hsc)), as.numeric(tclvalue(hsc)), as.numeric(tclvalue(hsc))-0.1), vscale = switch(ws, as.numeric(tclvalue(vsc))-0.08, as.numeric(tclvalue(vsc))-0.4, as.numeric(tclvalue(vsc))-0.4)))
	
	tkconfigure(ts_MZ_begin, state = "disable")
	tkconfigure(ts_MZ_range, state = "disable")
	
	tkconfigure(ts_v_angle, state = "disable")
	tkconfigure(ts_h_angle, state = "disable")
	
	tmpcl <- as.list(cl)
	tmpl_rt_begin <- list(substitute(as.numeric(tclvalue(vname)), list(vname = as.character(RT_begin))))
	tmpl_rt_range <- list(substitute(as.numeric(tclvalue(vname)), list(vname = as.character(RT_range))))
	
	
	
	

	
	
	names(tmpl_rt_begin) <- "rt_begin"
	names(tmpl_rt_range) <- "rt_range"
	
	
	cl <<- as.call(c(tmpcl, tmpl_rt_begin, tmpl_rt_range))
	exargs <<- c(exargs, tmpl_rt_begin, tmpl_rt_range)

    tkrreplot(img, hscale = switch(ws, as.numeric(tclvalue(hsc)), as.numeric(tclvalue(hsc)), as.numeric(tclvalue(hsc))-0.1), vscale = switch(ws, as.numeric(tclvalue(vsc))-0.08, as.numeric(tclvalue(vsc))-0.4, as.numeric(tclvalue(vsc))-0.4))
	tkconfigure(refresh.but, command = function() {
		tkrreplot(img,hscale = switch(ws, as.numeric(tclvalue(hsc)), as.numeric(tclvalue(hsc)), as.numeric(tclvalue(hsc))-0.1), vscale = switch(ws, as.numeric(tclvalue(vsc))-0.08, as.numeric(tclvalue(vsc))-0.4, as.numeric(tclvalue(vsc))-0.4))
	})
	if(length(NAME)==1) {tkbind(img, "<Double-Button-1>",OnLeftClick);tkconfigure(img,cursor="hand2")}
	tkconfigure(tt,cursor="arrow")
}





mzxml_tree <- function(outpath_temp){
	rootNode = paste(gsub(".*/(.*?)", "\\1",outpath_temp), "Node", sep = "")
	if(!inherits(try(tkinsert(treeWidget, "end", "root", rootNode, text = gsub(".*/(.*?)", "\\1",outpath_temp)), TRUE), "try-error")){
		ind_file = dir(outpath_temp, full.names = T, pattern = ".mzXML")
		ind_name = dir(outpath_temp, pattern = ".mzXML")
		ind_rep_name = unlist(lapply(ind_name, function(x) gsub("(.*?).mzXML", "\\1", x)))
		ind_name = unlist(lapply(ind_name, function(x) gsub("(.*?)_\\d.mzXML","\\1",x)))
		rep_time = table(ind_name)
		rep_time_cum = cumsum(rep_time)
		ind_name = unique(ind_name)
		for(ind in 1:length(ind_name)){
			name <- ind_name[ind]
			raw_list = if(ind==1) 1:rep_time_cum[ind] else (rep_time_cum[ind-1]+1):rep_time_cum[ind]

			nameNode <- ifelse(rep_time[ind]==1,paste(ind_file[raw_list],",,",sep=""),paste(ind_file[raw_list], collapse = ",,", sep = ""))
			tkinsert(treeWidget, "end", rootNode, nameNode, text = name)
			for(ind_rep in raw_list){
				repNode = ind_file[ind_rep]
				tkinsert(treeWidget, "end", nameNode, repNode, text = ind_rep_name[ind_rep])
			}
		}
		editPopupMenu <- tkmenu(treeWidget, tearoff = FALSE)
		tkadd(editPopupMenu, "command", label="2D Plot", command = function(){
			viewer <<- "2D"
			tkconfigure(save.but, state = "normal")
			eval(show2D())
		})
		tkadd(editPopupMenu, "command", label="3D Plot", command = function(){
			viewer <<- "3D"
			tkconfigure(save.but, state = "normal")
			eval(show3D())
		})
		tkadd(editPopupMenu, "command", label="TIC Plot", command = function(){
			viewer <<- "TIC"
			tkconfigure(save.but, state = "normal")
			eval(showTIC())
		})
		
		RightClick <- function(x, y) { 
			if(length(unlist(strsplit(tclvalue(tcl(treeWidget,"selection","get")), split = " |,,")))>1){
				tkentryconfigure(editPopupMenu, 0, state = "disable")
				tkentryconfigure(editPopupMenu, 1, state = "disable")
			}else if(length(unlist(strsplit(tclvalue(tcl(treeWidget,"selection","get")), split = " |,,")))==1){
				tkentryconfigure(editPopupMenu, 0, state = "normal")
				tkentryconfigure(editPopupMenu, 1, state = "normal")
			}else{
				return()
			}
			rootx <- as.integer(tkwinfo("rootx", treeWidget))  
			rooty <- as.integer(tkwinfo("rooty", treeWidget))
			xTxt <- as.integer(x) + rootx
			yTxt <- as.integer(y) + rooty
			
			.Tcl(paste("tk_popup", .Tcl.args(editPopupMenu, xTxt, yTxt)))
		}
		tcl(treeWidget, "bindArea", "<Button-3>", RightClick)
	}else{
		tkmessageBox(title = "Error", message = "Names are overlapping.\nPlease change the folder name, the folder names overlap.", icon = "error", type = "ok")
		cat("Data Import(Error) - Names are overlapping. Please change the folder name, the folder names overlap.\n")
	}
} 


DataImport <- function()
{
	dlg <- tktoplevel(); if(isIcon) tk2ico.set(dlg,icon)
	tkwm.title(dlg, "Data Import")

	fr_datatype <- tkframe(dlg)
	datatype <- tclVar("0")
	type_raw <- tkradiobutton(fr_datatype, variable = datatype, value = "1", command = function(){
		
		
		
		
		
		tkconfigure(inputlabel, text = "   Input folder (Raw data):      ")
		tkconfigure(box.input, command = function() tclvalue(textdatainput) <- tkchooseDirectory(initialdir = as.character(tclvalue(textdatainput))))
		tkconfigure(OK.but, text = "  Import  ")
		tkgrid.remove(fr_align3)
		tkgrid(fr_align1, sticky = "w")
		tkgrid(ProteoWizard_lab,ProteoWizard_widget,box.ProteoWizard)
		tkconfigure(ProteoWizard_lab,text="   Path of ProteoWizard:         ")
	})
	type_mzxml <- tkradiobutton(fr_datatype, variable = datatype, value = "2", command = function(){
		tkconfigure(box.input, state = "normal")
		tkconfigure(box.output, state = "normal")
		tkconfigure(textinputWidget, state = "normal")
		tkconfigure(textoutputWidget, state = "normal")
		tkconfigure(OK.but, state = "normal")
		tkconfigure(inputlabel, text = "   Input folder (mzXML file):      ")
		tkconfigure(box.input, command = function() tclvalue(textdatainput) <- tkchooseDirectory(initialdir = as.character(tclvalue(textdatainput))))
		tkconfigure(OK.but, text = "  Import  ")
		tkgrid.remove(fr_align1)
		tkgrid.remove(fr_align3)
	})
	type_align <- tkradiobutton(fr_datatype, variable = datatype, value = "3", command = function(){
		tkconfigure(box.input, state = "normal")
		tkconfigure(box.output, state = "normal")
		tkconfigure(textinputWidget, state = "normal")
		tkconfigure(textoutputWidget, state = "normal")
		tkconfigure(OK.but, state = "normal")
		tkconfigure(inputlabel, text = "   Input file:           ")
		tkconfigure(box.input, command = function() tclvalue(textdatainput) <- tkgetOpenFile(initialfile = as.character(tclvalue(textdatainput)), filetypes = "{{Text Files} {.txt .csv}}"))
		tkconfigure(OK.but, text = "  Continue  ")
		tkgrid.remove(fr_align1)
		tkgrid(fr_align3, sticky = "w")
		tkgrid(missing_lab, missing_widget, sticky = "w")
	})
	tkgrid(tklabel(fr_datatype,text="", font = fontIntro_para, height = 0))
	tkgrid(tklabel(fr_datatype, text = "   Data type:          "), type_raw, tklabel(fr_datatype, text = "Raw format", state=ifelse(osIsWin,"normal","disabled")), type_mzxml, tklabel(fr_datatype, text = "mzXML format"), type_align, tklabel(fr_datatype, text = "Peak abundance data   "), sticky = "w")
	tkgrid(fr_datatype, sticky = "w")
	tkconfigure(type_raw,state=ifelse(osIsWin,"normal","disabled"))
	fr_input <- tkframe(dlg)
	if(!exists("textdatainput")) textdatainput <<- tclVar("")
	textinputWidget <- tkentry(fr_input,width="55", textvariable = textdatainput, bg = "white")
	tkconfigure(textinputWidget, state = "disable")
	box.input <- tkbutton(fr_input, text = "...")
	tkconfigure(box.input, state = "disable")
	
	inputlabel <- tklabel(fr_input,text="   Input folder:        ")
	tkgrid(inputlabel, textinputWidget, box.input, tklabel(fr_input,text="    "), sticky = "w")
	
	if(!exists("textoutput")) textoutput <<- tclVar("")
	textoutputWidget <- tkentry(fr_input,width="55", textvariable = textoutput, bg = "white")
	tkconfigure(textoutputWidget, state = "disable")
	box.output <- tkbutton(fr_input, text = "...",  command = function() tclvalue(textoutput) <- tkchooseDirectory(initialdir = as.character(tclvalue(textoutput))))
	tkconfigure(box.output, state = "disable")
	tkgrid(tklabel(fr_input,text="   Output folder:    "), textoutputWidget, box.output, tklabel(fr_input,text="    "), sticky = "w")

	tkgrid(fr_input)
	
	fr_align1 <- tkframe(dlg)
	tkgrid(fr_align1, sticky = "w")
	text_ProteoWizard <<- tclVar("C:/Program Files (x86)/ProteoWizard")
	ProteoWizard_lab <- tklabel(fr_align1, text = "   Path of ProteoWizard:        ")
	ProteoWizard_widget <- tkentry(fr_align1, width = "55", textvariable = text_ProteoWizard, bg = "white", xscrollcommand=function(...){
	 msconvert <<- dir(tclvalue(text_ProteoWizard), full.names = T)
	 msconvert <<- msconvert[length(msconvert)]

	 if(file.exists(paste(msconvert,"/msconvert.exe",sep="/"))){
	  if(Sys.info()["sysname"]=="Windows"){
	   cat("ProteoWizard:",gsub(".*/","",msconvert),unique(ifelse(grepl("x86",msconvert),"32-bit",ifelse(grepl("x64",Sys.info()["release"]),"64-bit","32-bit"))),"\n")
	  }

	  tkconfigure(box.input, state = "normal")
	  tkconfigure(box.output, state = "normal")
	  tkconfigure(textinputWidget, state = "normal")
	  tkconfigure(textoutputWidget, state = "normal")
	  tkconfigure(OK.but, state = "normal")
	 } else {
	  tkmessageBox(title = "Error", message = "There is no 'msconvert.exe' under the specified folder.\nPlease check the path of ProteoWizard or check if ProteoWizard is installed successfully.", icon = "error", type = "ok")
	  cat("Data Import(Error) - There is no 'msconvert.exe' under the specified folder. Please check the path of ProteoWizard or check if ProteoWizard is installed successfully.\n")

	  tkconfigure(box.input, state = "disabled")
	  tkconfigure(box.output, state = "disabled")
	  tkconfigure(textinputWidget, state = "disabled")
	  tkconfigure(textoutputWidget, state = "disabled")
	  tkconfigure(OK.but, state = "disabled")
	 }
	})
	box.ProteoWizard <- tkbutton(fr_align1, text = "...",  command = function(){
	 tclvalue(text_ProteoWizard) <- tkchooseDirectory(initialdir = as.character(tclvalue(text_ProteoWizard)))
	})
	

	fr_align3 <- tkframe(dlg)
	tkgrid(fr_align3, sticky = "w")
	text_missing <<- tclVar("NA")
	missing_lab <- tklabel(fr_align3, text = "   Missing data:       ")
	missing_widget <- tkentry(fr_align3, width = 7, textvariable = text_missing, bg = "white")
	

	
	
	onImport_raw <- function()
	{
		if(file.exists(tclvalue(textdatainput)) & file.exists(tclvalue(textoutput))){
			tkdestroy(dlg)
			tkconfigure(tt,cursor="watch")
			pb <- tkProgressBar("Data Import - Please wait for data importing", "0% done", 0, 100, 0, width = 500)
			cat("Data Import - Please wait for data importing 0 ")
			inpath <- as.character(tclvalue(textdatainput))
			
			
			outpath <- tclvalue(textoutput)
			
			
			
			
			setwd(outpath)
			data_file = dir(inpath, full.names = T)
			data_name = dir(inpath)
			data_name = gsub("(.*?).raw","\\1",data_name)
			mzXML_path <- paste(outpath, "/mzXML", sep = "")
			dir.create(mzXML_path, showWarnings = F)
			textmzdatainput <<- tclVar(mzXML_path)
			
			
			for(x in 1:length(data_file)){
				system(paste('\"', msconvert, '/msconvert.exe\"', ' msconvert ', data_file[x], ' --32 --mzXML --filter "sortByScanTime"', ' -o ', mzXML_path, sep = ""), show.output.on.console = F)
				
				
				
				
				
				info <- sprintf("%d%%", round(100*x/length(data_file)))
				setTkProgressBar(pb, value = round(100*x/length(data_file)), sprintf("Data Import - Please wait for data importing (%s done)", info), paste(data_name[x], info, sep = " "))
				if(round(100*x/length(data_file))<100){
					cat(round(100*x/length(data_file)), " ", sep = "")
				}else{
					cat(round(100*x/length(data_file)), " \n", sep = "")
				}
			}
			
			
			
			setTkProgressBar(pb, value = 100, "Data Import - Please wait for data importing (100% done)", "Finished 100%")
			mzxml_tree(mzXML_path)
			Sys.sleep(1)
			close(pb)
			tkfocus(tt)
			tkentryconfigure(viewMenu, 0, state = "active")
			tkentryconfigure(viewMenu, 1, state = "active")
			tkentryconfigure(viewMenu, 2, state = "active")
			tkentryconfigure(viewMenu, 3, state = "active")
			tkentryconfigure(peakMenu, 0, state = "active")
			tkconfigure(tt,cursor="arrow")
			tkmessageBox(title="Data Import",message="Data import is done.", icon = "info", type = "ok")
			cat("Data Import-Data import is done.\n", sep = "")
		}else{
			tkmessageBox(title = "Error", message = "File is not found.\nPlease input the correct directories.", icon = "error", type = "ok")
			cat("Data Import(Error) - File is not found. Please input the correct directories.\n")
			tkfocus(dlg)
		}
	}
	
	onImport_mzxml <- function()
	{
		if(file.exists(tclvalue(textdatainput)) & file.exists(tclvalue(textoutput))){
			if(sum(grepl(".mzXML || .d", dir(tclvalue(textdatainput))))){
				tkdestroy(dlg)
				tkconfigure(tt,cursor="watch")
				
				inpath <- as.character(tclvalue(textdatainput))
				setwd(tclvalue(textoutput))
				textmzdatainput <<- tclVar(inpath)
				mzxml_tree(inpath)
				
				
				tkfocus(tt)
				tkentryconfigure(viewMenu, 0, state = "active")
				tkentryconfigure(viewMenu, 1, state = "active")
				tkentryconfigure(viewMenu, 2, state = "active")
				tkentryconfigure(viewMenu, 3, state = "active")
				tkentryconfigure(peakMenu, 0, state = "active")
				tkconfigure(tt,cursor="arrow")
				tkmessageBox(title="Data Import",message="Data import is done.")
				cat("Data Import-Data import is done.\n", sep = "")
			}else{
				tkmessageBox(title = "Error", message = "Wrong file format.\nPlease check the type of all files in the input folder is in mzXML format.", icon = "error", type = "ok")
				cat("Data Import(Error)-Wrong file format. Please check the type of all files in the input folder is in mzXML format.\n")
				tkfocus(dlg)
			}
		}else{
			tkmessageBox(title = "Error", message = "File is not found.\nPlease input the correct directories.", icon = "error", type = "ok")
			cat("Data Import(Error) - File is not found. Please input the correct directories.\n")
			tkfocus(dlg)
		}
	}
	
	onImport_align <- function()
	{
		if(file.exists(tclvalue(textdatainput)) & file.exists(tclvalue(textoutput))){
			tkdestroy(dlg)
			checking <- tktoplevel(); if(isIcon) tk2ico.set(checking,icon)
			tktitle(checking) = "Column Information"
			input_idx_head = tkframe(checking)
			tkgrid(tklabel(input_idx_head, text = "   Please provide column information:"), sticky = "w")
			tkgrid(tklabel(input_idx_head, text = "", font = fontIntro_para, height = 0), sticky = "w")
			tkgrid(input_idx_head, sticky = "w")
			
			input_idx = tkframe(checking)
			
			PI.lab = tklabel(input_idx, text = "   Column index of peak index: ")
			tk2tip(PI.lab, "Please input column index of peak index if you have.")
			index.PI <<- tclVar("NA")
			PI_entry = tkentry(input_idx, width = "5", textvariable = index.PI, bg = "white")	
			tkgrid(PI.lab, PI_entry, sticky = "w")
			
			
			mass.lab <- tklabel(input_idx, text = "   Column index of mass: ")
			index.mass <<- tclVar("")
			mass_entry = tkentry(input_idx, width = "5", textvariable = index.mass, bg = "white")
			tkgrid(mass.lab, mass_entry, sticky = "w")
			
			
			RT.lab = tklabel(input_idx, text = "   Column index of retention time: ")
			index.RT <<- tclVar("")
			RT_entry = tkentry(input_idx, width = "5", textvariable = index.RT, bg = "white")
			RT_unit = c("min", "sec")
			RT_text = tclVar("sec")
			RTunit.lab = tklabel(input_idx, text = "  unit: ")
			RTcomboBox <- ttkcombobox(input_idx, state = "readonly", values = RT_unit, width = 4, textvariable = RT_text)
			tkgrid(RT.lab, RT_entry, RTunit.lab, RTcomboBox, tklabel(input_idx, text = "   "), sticky = "w")
			
			
			first.lab = tklabel(input_idx, text = "   Column index of first sample: ")
			index.first <<- tclVar("")
			first_entry = tkentry(input_idx, width = "5" ,textvariable = index.first, bg = "white")
			tkgrid(first.lab, first_entry, sticky = "w")
			
			tkgrid(tklabel(input_idx, text = "", height = 0, font = fontIntro_para))
			tkgrid(input_idx)
			
			contiframe = tkframe(checking)

			conti.button <- tkbutton(contiframe, text = "   Import   ", command = function(){
				RT.unit <<- tclvalue(RT_text)
				tkdestroy(checking)
				tkconfigure(tt,cursor="watch")
				
				align_table <- read.table(tclvalue(textdatainput), header = T, quote = "\"", fill = T, sep = ifelse(grepl(".txt", tclvalue(textdatainput)), "\t", ","), as.is = T, na.string = c("NA", as.character(tclvalue(text_missing))))
				if(as.character(tclvalue(index.PI))=="NA"){
					align_table = cbind(1:nrow(align_table), align_table[,as.numeric(c(tclvalue(index.mass), tclvalue(index.RT), seq.int(as.numeric(tclvalue(index.first)), ncol(align_table))))])
				}else{
					align_table = align_table[,as.numeric(c(tclvalue(index.PI), tclvalue(index.mass), tclvalue(index.RT), seq.int(as.numeric(tclvalue(index.first)), ncol(align_table))))]
				}
				nRow <- nrow(align_table)
				nCol <- ncol(align_table)
				align_table <- peakabunCheck(align_table)
					
					
				if(RT.unit == "min") align_table[,3] = 60*align_table[,3]
				colnames(align_table) = c("Peak_Index", "mz", "Ret_time.sec", gsub("^\\d{4,8}_(.*?)|(^[A-Z]\\d{4,8})_(.*?)", "\\1", colnames(align_table)[-c(1:3)]))
				write.table(align_table, file = gsub(".csv|.txt", "_forSMART.csv", tclvalue(textdatainput)), sep = ",", col.names = T, row.names = F, quote = c(1))
				textAbuninput <<- tclVar(gsub(".csv|.txt", "_forSMART.csv", tclvalue(textdatainput)))
				setwd(tclvalue(textoutput))
				tkfocus(tt)
				
				tkentryconfigure(peakMenu, 1, state = "active")
				tkentryconfigure(QCsubMenu, 0, state = "active")
				
				tkentryconfigure(batchsubMenu, 0, state = "active")
				tkentryconfigure(batchsubMenu, 1, state = "active")
				tkentryconfigure(statMenu, 1, state = "active")
				tkentryconfigure(statMenu, 2, state = "active")
				tkconfigure(tt,cursor="arrow")

				tkmessageBox(title = "Data Import", message = sprintf("Data import is done.\n%d out of %d non-informative replicate sample(s) are excluded.\n%d out of %d non-informative peak(s) are excluded.",nCol-ncol(align_table),nCol-3,nRow-nrow(align_table),nRow), icon = "info", type = "ok")
				cat(sprintf("Data import is done.\n%d out of %d non-informative replicate sample(s) are excluded.\n%d out of %d non-informative peak(s) are excluded.\n",nCol-ncol(align_table),nCol-3,nRow-nrow(align_table),nRow))
				
				
			})
			tkgrid(conti.button)
			tkgrid(tklabel(contiframe, text = "", height = 0, font = fontIntro_para))
			tkgrid(contiframe)
			
		}else{
			tkmessageBox(title = "Error", message = "File is not found.\nPlease input the correct directories.", icon = "error", type = "ok")
			cat("Data Import(Error) - File is not found. Please input the correct directories.\n")
			tkfocus(dlg)
		}
	}
	onImport <- function(){
		cat(sprintf("Data Import - Input: %s.\n",tclvalue(textdatainput)))
		cat(sprintf("Data Import - Output: %s.\n",tclvalue(textoutput)))
		if(tclvalue(datatype)=="1"){
			eval(onImport_raw())
		}else if(tclvalue(datatype)=="2"){
			eval(onImport_mzxml())
		}else{
			eval(onImport_align())
		}
	}
	onCancel <- function()
	{
		
		tkdestroy(dlg)
		tkfocus(tt)
	}
	fr <- tkframe(dlg)
	tkgrid(tklabel(fr,text="", font = fontIntro_para, height = 0))
	OK.but     <-tkbutton(fr,text="  Import  ",command=onImport)
	tkconfigure(OK.but, state = "disable")
	Cancel.but <-tkbutton(fr,text="  Cancel  ",command=onCancel)
	tkgrid(tklabel(fr,text="   "), OK.but,tklabel(fr,text="                                                         "), Cancel.but, tklabel(fr,text="   "))
	tkgrid(tklabel(fr,text="    ", font = fontIntro_para, height = 0))
	tkgrid(fr)
	

	tkfocus(dlg)
	tkbind(dlg, "<Destroy>", function() {tkgrab.release(dlg);tkfocus(tt)})
	
	
	tkwait.window(dlg)

	
} 

r.value = function(temp.data, num.replicate, index.group) {
	result = t(sapply(1:nrow(temp.data), function(i) {
		data.x = unlist(temp.data[i,])
		
		y = tapply(data.x, index.group, function(x) if(length(x) == num.replicate) x else c(x, rep(NA, num.replicate-length(x)))) 
		data.matrix = do.call(cbind, y)
		
		temp = colMeans(data.matrix, na.rm = T)
		
		individual = mean(apply(data.matrix, 2, var, na.rm = T), na.rm = T)
		individual.cv = mean(apply(data.matrix, 2, sd, na.rm = T)/temp, na.rm = T)
		
		group = var(temp, na.rm = T)
		group.cv = sd(temp, na.rm = T)/mean(temp, na.rm = T)
		
		c(Within = individual, Between = group, Within_cv = individual.cv, Between_cv = group.cv)
	}))
	return(ifelse(result[,2]==0, ifelse(result[,1]==0, 0, 10000), result[,1]/result[,2]))
}  

r.cluster = function(Phx.E, r = r, cut.r = Inf, num.replicate, sample.info.temp, c.path, index.plot = T, index.plot.sh = T, index.correlation = T, index.output = F, metric = "pearson", method, height_method, QC = T){
	if(index.correlation){
		data.temp = Phx.E[(r < cut.r),]
		sd.temp = apply(data.temp, 2, sd, na.rm=T)
		index.sd.0 = which(sd.temp == 0)
		if(length(index.sd.0)){
			data.temp[, index.sd.0] = data.temp[,index.sd.0] + runif(length(index.sd.0)*nrow(data.temp), 0, 10^-5)
		}
		cor.Phx = cor(data.temp, method = metric, use = "p")
		rownames(cor.Phx) = colnames(cor.Phx) = sample.info.temp$Group 
	} else {
		cor.Phx = Phx.E[(r < cut.r),]
		colnames(cor.Phx) = sample.info.temp$Sample 
	}
	subtitle = ifelse(QC, paste(sum(r < cut.r), " peaks with r-value < ", cut.r, sep = ""), paste(sum(r < cut.r), " Scans", sep = ""))
	
	bcc = clust.all.2(cor.Phx, metric = metric, method = method, index.correlation = index.correlation, operation = paste("r_", cut.r, sep = ""), subtitle = subtitle, path = c.path, num.replicate = num.replicate, num.sample = length(unique(sample.info.temp$Group)), index.plot = index.plot, index.plot.sh = index.plot.sh, sample.barcode = sample.info.temp$Sample, sample.id = sample.info.temp$Group, index.output = index.output, height_method = height_method, QC = QC)

	bcc$num.probe = sum((r < cut.r))
	ifelse(index.output, return(bcc), return())
}

my.cluster = function(x, m) { 
	x = -x
	temp = m[x, ]
	temp = -temp
	for(i in temp){
		if(i > 0){
			result <<- c(result, i)
		} else Recall(i, m)
	}
}

clust.all.2 = function(x, index.correlation = T, metric = "pearson", method = c("complete", "average", "ward"), operation = "ITRI", index.plot = T, num.replicate, sample.barcode = NULL, subtitle = "All Peaks", path, num.sample = NULL, sample.id = NULL, d.num = 120, index.output, height_method, index.plot.sh, QC){
	
	
	

	bcc <- vector("list",14); names(bcc) <- c("accuracy","height","height_thresh","steps","point","backcol","index.bad","index.bad_label","index.h2","index.h2_label","tree","col.cluster","sb.bad","sb.bad2")
	for(i in names(bcc)){
		if(i=="accuracy") bcc[[i]] <- as.vector(rep(NA,length(method)),"list")
		else bcc[[i]] <- vector("list",length(method))
		names(bcc[[i]]) <- method
	}

	if(index.correlation) xx <- as.dist((1-x)/2) 
	else xx <- dist(t(x)) 

	if(any(is.na(xx))) return(bcc) 
	
	for(m in method){
		clust.temp <- hclust(xx, ifelse(m == "ward", paste(m, ".D", sep = ""), m)) 
		
		title.temp = paste(toupper(substring(m, 1,1)), substring(m, 2), sep = "")
		if(QC){
			path.temp <- sprintf("%s/%s/%s/%s", path, c("Parametric", "Nonparametric", "Bootstrap")[height_method], paste0(toupper(substring(metric, 1, 1)), substring(metric, 2)), title.temp)
			dir.create(path.temp, recursive = T, showWarnings = F)
			
			
			
			dir.create(path.temp, showWarnings = F)
		}else{
			path.temp = paste(path, "/", title.temp, sep = "")
			dir.create(path.temp, showWarnings = F)
		}

		
		if(m == "ward") {
			title = paste(title.temp, "Method", sep = " ")
		} else title = paste(title.temp, "Linkage", sep = " ")
		
		
		merge.m = clust.temp$merge
		label.all = clust.temp$labels
		height.all = clust.temp$height
		
		clust.check = lapply(clust.temp$order, function(x){
			
			now = label.all[x] 
			
			x = -x
			n.row = (which(merge.m == x) - 1)%%(nrow(merge.m)) + 1 
			
			n.col = (which(merge.m == x) - 1)%/%(nrow(merge.m)) + 1 
			
			temp = -merge.m[n.row, ifelse(n.col == 1, 2, 1)]
			
			result <<- NULL
			if(temp > 0){ 
				result <<- c(result, temp)
			} else my.cluster(temp, merge.m)
			
			
			
			
			
			list(clust.result = result, index = as.integer(all(label.all[result] == now)), height = clust.temp$height[n.row], 
				row = n.row, column = n.col, order = -x)
		})
		names(clust.check) = label.all[clust.temp$order]
		

		
		clust.tree = sapply(clust.check, function(x) x[[1]])
		clust.height = sapply(clust.check, function(x) x[[3]])
		clust.row = sapply(clust.check, function(x) x[[4]])
		clust.col = sapply(clust.check, function(x) x[[5]])
		clust.order = sapply(clust.check, function(x) x[[6]])
	

		
		name.all = unique(label.all[clust.temp$order])
		determine.lowest = function(index, now, h, r, cc, m, o){
			ii = order(h)
			index = index[ii] 
			h = h[ii] 
			r = r[ii] 
			cc = cc[ii] 
			lowest = index[1] 
			for(i in 1:length(index)){
				n.row = r[i] 
				n.col = cc[i]
				temp = -m[n.row, ifelse(n.col == 1, 2, 1)] 
				
				if(temp > 0){ 
					
					if(temp %in% o){ 
						lowest = index[i]
						
						
						
						break
					}
				}
			}
			return(lowest)
		}
	
		
		test = function(index, index2, now, clust.row, clust.col, clust.order, merge.m, o.temp, label.all, height.all, sample.barcode){
			
			
		
			n.index = length(index)
			n.row = clust.row[index2] 
			
			n.col = clust.col[index2]
			
			
			
			
			
			
			index.order = clust.order[index2] 
			height = 0
			
			clust = T 
			step = 0 
			i = 1 
			s = 1 
			
			
			
			while(i != n.index){
				
				temp = -merge.m[n.row, ifelse(n.col == 1, 2, 1)] 
				
				result <<- NULL 
				
				
				if(temp > 0){ 
					
					if(temp %in% o.temp){ 
						i = i + 1 
						
						index.order = c(index.order, temp) 
						height = c(height, height.all[n.row]) 
						
						clust = c(clust, if((tail(step, 1) + 1) == s) tail(clust, 1) else F)
						step = c(step, s) 
					}
				} else { 
					my.cluster(temp, merge.m) 
					
					
					nn = sum(result %in% o.temp) 
					i = i + nn
					if(nn){
						
						index.order = c(index.order, result[result %in% o.temp])
						height = c(height, rep(height.all[n.row], nn)) 
						
						clust = c(clust, if((tail(step, 1) + 1) == s & length(result) == nn) rep(tail(clust, 1), nn) else rep(F, nn))
						step = c(step, rep(s, nn)) 
					}
					

					
				}
				
				
				s = s + 1
				
				

				k = n.row
				
				
				
				n.row = (which(merge.m == k) - 1)%%(nrow(merge.m)) + 1
				n.col = (which(merge.m == k) - 1)%/%(nrow(merge.m)) + 1
			}
			
			
			a = cbind(Order = index.order, OrderTree = index[match(index.order, o.temp)], Height = height, Step = step, Clust = clust)
			
			rownames(a) = sample.barcode[index.order]
			return(a)
		}
	
		clust.stack2 = lapply(name.all, function(x){
		
			index = which(names(clust.tree) %in% x) 
			l.temp = clust.tree[index] 
			h.temp = clust.height[index] 
			r.temp = clust.row[index] 
			c.temp = clust.col[index] 
			o.temp = clust.order[index] 
			
				
			index2 = determine.lowest(index, now = x, h = h.temp, r = r.temp, cc = c.temp, m = merge.m, o = o.temp) 
			
			
			
			a = test(index, index2, now = x, clust.row, clust.col, clust.order, merge.m, o.temp, 
				label.all, height.all, sample.barcode)
			
			return(a)
		})
		names(clust.stack2) = name.all
		index.major = sapply(clust.stack2, function(x) sum(x[, "Clust"] == 0) > floor(num.replicate/2)) 
		
		if(sum(index.major)) {
			clust.stack3 = lapply(which(index.major), function(i){
				a = clust.stack2[[i]]
				index = subset(a, subset = a[, "Clust"] == 0, select = OrderTree) 
				l.temp = clust.tree[index] 
				h.temp = clust.height[index] 
				r.temp = clust.row[index] 
				c.temp = clust.col[index] 
				o.temp = clust.order[index] 
				x = names(clust.stack2)[i]
				index2 = determine.lowest(index, now = x, h = h.temp, r = r.temp, cc = c.temp, m = merge.m, o = o.temp) 
				a = test(index, index2, now = x, clust.row, clust.col, clust.order, merge.m, o.temp,
					label.all, height.all, sample.barcode)
			})
			
			index.renew = sapply(clust.stack3, function(x) all(x[, "Clust"] == 1)) 
			
			if(sum(index.renew)){
				clust.stack4 = lapply(which(index.major)[index.renew], function(i){
					a = clust.stack2[[i]]
					
					index = subset(a, subset = a[, "Clust"] == 0, select = OrderTree) 
					h.temp = clust.height[index] 
					r.temp = clust.row[index] 
					c.temp = clust.col[index] 
					o.temp = clust.order[index]
					x = names(clust.stack2)[i]
					index2 = determine.lowest(index, now = x, h = h.temp, r = r.temp, cc = c.temp, m = merge.m, o = o.temp) 
					index = subset(a, select = OrderTree) 
					o.temp = clust.order[index] 
					a = test(index, index2, now = x, clust.row, clust.col, clust.order, merge.m, o.temp,
						label.all, height.all, sample.barcode)
					return(a)
				})
			
				
				for(i in 1:sum(index.renew)) clust.stack2[[which(index.major)[index.renew][i]]] = clust.stack4[[i]]
			}
		}
		sample_ID_rep = unlist(lapply(clust.stack2, rownames))
		s = unlist(lapply(clust.stack2, function(x) max(x[-1, "Step"]))) 
		h = unlist(lapply(clust.stack2, function(x) max(x[-1, "Height"]))) 
		h_full = unlist(lapply(clust.stack2, function(x){
			temp = x[-1, "Height"]
			c(temp[1], temp)
		}))
		s_full = unlist(lapply(clust.stack2, function(x){
			temp = x[-1, "Step"]
			c(temp[1], temp)
		}))
		index.good = sapply(clust.stack2, function(x) all(x[, "Clust"] == 1))
		clust.bad = unlist(lapply(clust.stack2, function(x) any(x[, "Clust"]==0))) 
		clust.col = unlist(lapply(clust.stack2, function(x){
			if(sum(x[, "Clust"])<(num.replicate*2/3)){
				return(2)
			}else if(sum(x[, "Clust"])>=(num.replicate*2/3) & sum(x[, "Clust"])<num.replicate){
				ifelse(all(x[,"Clust"]==1), return(1), return(3))
			}else{
				return(1)
			}
		}))
		clust.col_full = unlist(lapply(clust.stack2, function(x){
			if(sum(x[, "Clust"])<(num.replicate*2/3)){
				return(rep("poor", nrow(x)))
			}else if(sum(x[, "Clust"])>=(num.replicate*2/3) & sum(x[, "Clust"])<num.replicate){
				ifelse(all(x[,"Clust"]==1), return(rep("perfect", nrow(x))), return(rep("well", nrow(x))))
			}else{
				return(rep("perfect", nrow(x)))
			}
		}))
		index.bad = which(unlist(lapply(clust.stack2, function(x) if(any(x[, "Clust"]==0)){return(1)}else{return(0)}))==1)
		
			
			
			
		
		index.bad_label = unlist(lapply(clust.stack2, function(x){
			if(length(which(x[-1, "Clust"]==0))>0){
				if(sum(x[, "Clust"])<(num.replicate*2/3)){
					return(gsub("_[0-9]", "", names(which(x[, "Clust"]==0))[1]))
				}else{
					ind_split = unlist(strsplit(names(which(x[, "Clust"]==0)), split = "_"))
					
					temp = matrix(c(paste(ind_split[1:(length(ind_split)-1)], collapse = "_"), ind_split[length(ind_split)]), ncol = 2, byrow = T)
					x = clust.height[x[which(x[, "Clust"]==0), "OrderTree"]]
					temp1 = temp[1,1]
					for(i in 1:nrow(temp)){
						temp1 = paste(temp1, "_", temp[i, 2], sep = "")
						if(x[i] != x[i+1] & i < nrow(temp)){
							temp1 = paste0(temp1, "|")
						}
					}
					return(temp1)
				}
			}else{ return()}
		}))
		
		
		
		
		
			
		
		h.good = unlist(lapply(which(index.good), function(i) clust.stack2[[i]][-1, "Height"]))
		h.good[h.good==0] = min(h.good[h.good!=0])/10
		h.good.max = as.vector(sapply(which(index.good), function(i) max(clust.stack2[[i]][-1, "Height"])))
		
		col.cluster = unlist(lapply(clust.stack2, function(x){
			if(sum(x[, "Clust"])<(num.replicate*2/3)){
				return(1)
			}else if(sum(x[, "Clust"])>=(num.replicate*2/3) & sum(x[, "Clust"])<num.replicate){
				ifelse(all(x[,"Clust"]==1), return(3), return(2))
			}else{
				return(3)
			}
		}))[label.all[clust.temp$order]]
		
		good.height = h.good.max
		if(index.plot.sh & length(h.good)>1){ 
			
			h1 = log(h.good)
			h2 = log(h.good.max)
			pv = round(ks.test(h1, "pnorm", mean(h1), sd(h1))$p.value, 3)
			
			
			
			
			
			
			
			
			if(height_method == 1){  
				bd = mean(h1) + 3*sd(h1); 
			}else if(height_method == 2){  
				bd = quantile(h1, probs = 0.75) + 1.5*diff(quantile(h1, probs = c(0.25, 0.75)))
			}else if(height_method == 3){  
				b = boot(h1, function(u,i) mean(u[i])+3*sd(u[i]), R = 5000)
				bd = mean(boot.ci(b, type = "bca")$bca[4:5])
			}
			
			index.height = which(h > exp(bd))
			
			plot.target = function(index.bad, s, h, p.label, ...){
				temp = subset(rbind(s, h), select = index.bad)
				colnames(temp) = p.label
				temp = t(subset(temp, select = order(temp[1,], temp[2,])))
				temp1 = unique(temp)
				rt = rownames(temp)
				rt1 = rownames(temp1)
				ll <<- 0
				dup.name = sapply(1:nrow(temp1), function(i){
					j = ll + 2 
					if(j <= nrow(temp)){
						name = rt1[i]
						while(all(temp[j,] == temp1[i,])){
							name = c(name, rt[j])
							j = j + 1
							if(j > nrow(temp)) break
						}
						ll <<- length(name) + ll
						nn = paste(name, collapse = " & ")
						
						
						
							 
						
						
						
						
						return(nn)
					} else return(tail(rt1, 1))
				})
				text(temp1[,1], temp1[,2], labels = dup.name, pos = 3, ...)
				return(dup.name)
			}
			pch = c(21, 24)[(h > exp(bd)) + 1]
			png(filename = paste(path.temp, "/Step_vs_h_", operation, ".png", sep = ""), width = 1920, height = 998)
			
			par(mgp = c(2.5, 0.6, 0), cex.lab = 2, cex.main = 2, cex.axis = 1.5, tck = -0.005, font.lab = 2)

			
				
			plot(s[order(clust.col)], h[order(clust.col)], xlim = c(0, max(s) * 1.1), ylim = c(0, max(c(h, exp(bd)))*1.1), pch = pch[order(clust.col)], bg = sort(clust.col), lwd = 1,
				xlab = "Steps", ylab = "", cex = 2)
			mtext("h",2,line=2.5,las=2,cex = 2)
			title(main = paste("Poorly-clustered: ", sum(clust.col==2), " samples, ", gsub("_", ": ", operation), sep = ""))
			abline(h = exp(bd), col = 3)
			text(0, exp(bd), label = round(exp(bd), 3), col = 3, font = 4, cex = 3)

		
			
			sb.bad = sb.bad2 = NA
			
			if(length(index.bad)) sb.bad = plot.target(index.bad, s, h, index.bad_label, col = "gray38", cex = 1.4, xpd = NA, offset = 0.8)
			
			
			index.h2 = index.height[!index.height %in% index.bad]
			index.h2_label = unlist(lapply(clust.stack2[names(index.h2)], function(x){
				temp = matrix(unlist(strsplit(rownames(x)[which(x[,"Height"]>exp(bd))], split = "_")), ncol = 2, byrow = T)
				temp1 = temp[1,1]
				for(i in 1:nrow(temp)){temp1 = paste(temp1, "_", temp[i, 2], sep = "")}
				return(temp1)
			}))
			if(length(index.h2)) sb.bad2 = plot.target(index.h2, s, h, index.h2_label, col = "gray38", cex = 1.4, xpd = NA, offset = 0.8)
			legend(x = "topleft", legend = c(paste("h < ", round(exp(bd),3), sep = ""), paste("h > ", round(exp(bd),3), sep = ""), "perfectly-clustered", "well-clustered", "poorly-clustered"), col = c(1, 1, 1, 3, 2), pch = c(21, 24, 16, 16, 16), bty = "n", cex = 2, inset = c(0,0))
			dev.off()
			
			
			
			
			
			
			
			sb = sample.barcode[clust.temp$order]
			cb = unlist(lapply(clust.stack2, function(x) x[, "Clust"])) 
			names(cb) = unlist(lapply(clust.stack2, function(x) rownames(x)))
			
			
		}

		if(index.plot){
			
			if(QC){
				png(paste(path.temp, "/Clustering_", operation, ".png", sep = ""), width = 1920, height = 1080)
			}else{
				png(paste(path.temp, "/Clustering.png", sep = ""), width = 1920, height = 1080)
			}
			
			d = d.num 
			num.separate = ceiling(length(col.cluster)/d)

			remainder = (length(col.cluster) - 1) %% d + 1
			ct = as.dendrogram(clust.temp)
			
			
				
					
					
				
				
			
			
			k = length(col.cluster)-1

			layout(matrix(c(1, (1:(2*num.separate)) + 1, 0), ncol = 1), heights = c(0.7, rep(c(1.5, 0.05), num.separate), 0.05))
			
			
			
			par(mar = c(0, 0, 0, 0))
			plot(1, type = "n", axes = F)
			
			title(paste("\n\n", title, " (", subtitle, ")", sep = ""), cex.main = 3);
			legend("bottomleft", legend = sprintf("%s (%d)",c("perfectly-clustered sample", "well-clustered sample", "poorly-clustered sample"),table(factor(clust.col,levels=1:3))[c("1","3","2")]), pch = 15, bty = "n", horiz = T, col = c("blue", "green", "red"), pt.cex = 2, cex = 2)

			for(i in 1:num.separate){
			
				
				par(mar = c(3.8, 0, 0, ifelse(i==num.separate, (d-remainder)/d*par("din")[1], 0)))
				ct1 = ct
				plot(ct, xaxs = "i", axes = F, nodePar=list(pch = c(NA, NA), cex = 0.8, lab.cex = 1.2, col = 3), 
					dLeaf = 0.04*max(clust.temp$height),
					
					xlim = 0.5 + c(d*(i-1), ifelse(i == num.separate, d*(i-1) + remainder, d*i)), edgePar = list(lwd = 2))

				

				par(mar = c(0, 0, 0, ifelse(i==num.separate, (d-remainder)/d*par("din")[1], 0)))
				

				image(cbind(1:length(col.cluster)), col = c("red", "green", "blue")[col.cluster], axes = FALSE,
					xlim = -1/k/2 + c(d*(i-1)/k, ifelse(i == num.separate, (d*(i-1) + remainder)/k, d*(i)/k))
						
				)

				
			}
			dev.off()
		} 
		
		if(index.output & length(h.good)>1){ 
			
			bcc$height[[m]] = h
			bcc$height_thresh[[m]] = exp(bd)
			bcc$steps[[m]] = s
			bcc$point[[m]] = pch
			bcc$backcol[[m]] = clust.col
			bcc$index.bad[[m]] = index.bad
			bcc$index.bad_label[[m]] = index.bad_label
			bcc$index.h2[[m]] = index.h2
			bcc$index.h2_label[[m]] = index.h2_label
			bcc$tree[[m]] = ct
			bcc$col.cluster[[m]] = col.cluster
			
			if(length(index.bad)) bcc$sb.bad[[m]] = sb.bad
			if(length(index.h2)) bcc$sb.bad2[[m]] = sb.bad2
			
			
			write.table("=====================", file = paste(path.temp, "/QC_", operation, ".csv", sep = ""), row.names = F, col.names = F, quote = F, sep = ",")
			write.table(paste("S-value threshold is ", round(exp(bd), 4), sep = ""), file = paste(path.temp, "/QC_", operation, ".csv", sep = ""), row.names = F, col.names = F, quote = F, sep = ",", append = T)
			write.table("=====================", file = paste(path.temp, "/QC_", operation, ".csv", sep = ""), row.names = F, col.names = F, quote = F, sep = ",", append = T)
			
			index.bad_delete = rep("No", length(sample_ID_rep))
			index.bad_delete[which((clust.col_full=="poor") | ((clust.col_full == "well")&(h_full>exp(bd))))] = "Yes"
			
			bcc$accuracy[[m]] = mean(ifelse(index.bad_delete=="Yes", 0, 1))
			
			result_output = data.frame(SampleID_rep = sample_ID_rep, remove_sample = index.bad_delete, Quality = clust.col_full, h = round(h_full, 4), stringsAsFactors = F)
			sorting = data.frame(h_full, 1:length(h_full), sample_ID_rep, as.numeric(factor(gsub("_[0-9]", "", sample_ID_rep))))
			poor_sort = sorting[which(clust.col_full=="poor"),]
			poor_sort_index = poor_sort[order(poor_sort[,1], decreasing = T), 2]
			well_yes_sort = which(sorting[,4] %in% unique(sorting[clust.col_full == "well" & index.bad_delete == "Yes",4]))
			well_yes_sort = sorting[well_yes_sort[clust.col_full[well_yes_sort] == "well"],]
			well_yes_sort_index = well_yes_sort[order(well_yes_sort[,1], decreasing = T), 2]
			well_no_sort = which(!sorting[,4] %in% unique(sorting[clust.col_full == "well" & index.bad_delete == "Yes",4]))
			well_no_sort = sorting[well_no_sort[clust.col_full[well_no_sort] == "well"],]
			well_no_sort_index = well_no_sort[mixedorder(well_no_sort[,3]), 2]
			perfect_sort = sorting[which(clust.col_full == "perfect"),]
			perfect_sort_index = perfect_sort[mixedorder(perfect_sort[,3]), 2]
			result_output = result_output[c(poor_sort_index, well_yes_sort_index, well_no_sort_index, perfect_sort_index),]
			
			
			write.table(result_output, file = paste(path.temp, "/QC_", operation, ".csv", sep = ""), row.names = F, col.names = T, quote = F, sep = ",", append = T)
			
			
			
			
			
			
			
		}
	}
	ifelse(index.output, return(bcc), return()) 
} 




TICclustering <- function()
{
	dlg <- tktoplevel(); if(isIcon) tk2ico.set(dlg,icon)
	tkwm.title(dlg, "TIC Clustering")
	fr_input <- tkframe(dlg)
	
	
		
	
	
	
	
	tkgrid(tklabel(fr_input, text = "", height = 0, font = fontIntro_para))
	
	if(file.exists(paste(tclvalue(textoutput), "/TIC Clustering/TICtable.csv", sep = ""))){
		tkgrid(tklabel(fr_input, text = paste("Input file: ", tclvalue(textoutput), "/TIC Clustering/TICtable.csv", "    ", sep = ""), padx = 11, wraplength = 600, justify = "left"), sticky = "w")
	}else{
		tkgrid(tklabel(fr_input, text = paste("Input folder: ", tclvalue(textmzdatainput), "    ", sep = ""), padx = 11, wraplength = 600, justify = "left"), sticky = "w")
	}
	
	
	
	tkgrid(tklabel(fr_input, text = paste("Output folder: ", tclvalue(textoutput), "    ", sep = ""), padx = 11, wraplength = 600, justify = "left"), sticky = "w")
	tkgrid(tklabel(fr_input, text = "", height = 0, font = fontIntro_para))
	tkgrid(fr_input, sticky = "w")
	
	fr_RT <- tkframe(dlg)
	RT_start <- tclVar("0")
	rt_start.entry <- tkentry(fr_RT, width = "7", textvariable = RT_start, bg = "white")
	RT_end <- tclVar("600")
	rt_end.entry <- tkentry(fr_RT, width = "7", textvariable = RT_end, bg = "white")
	tkgrid(tklabel(fr_RT, text = "   Retention time:      From"), rt_start.entry, tklabel(fr_RT, text = "      to"), rt_end.entry, tklabel(fr_RT, text = "   (sec)"), sticky = "w")
	tkgrid(fr_RT, sticky = "w")
	
	fr_method <- tkframe(dlg)
	method_c.val <- tclVar("0")
	method_c <- tkcheckbutton(fr_method, variable = method_c.val)
	method_a.val <- tclVar("1")
	method_a <- tkcheckbutton(fr_method, variable = method_a.val)
	method_w.val <- tclVar("0")
	method_w <- tkcheckbutton(fr_method, variable = method_w.val)
	tkgrid(tklabel(fr_method, text = "   Clustering method:     "), method_c, tklabel(fr_method, text = "Complete linkage"), method_a, tklabel(fr_method, text = "Average linkage"), method_w, tklabel(fr_method, text = "Ward's method                                                                                     "), sticky = "w")
	tkgrid(tklabel(fr_method, text = "", height = 0, font = fontIntro_para))
	tkgrid(fr_method, sticky = "w")
	
	
	
	
	
	
	
	
	
	
	
	onClust <- function()
	{
		tkdestroy(dlg)
		tkconfigure(tt,cursor="watch")
		inpath = tclvalue(textmzdatainput)
		outpath = paste(tclvalue(textoutput), "/TIC Clustering", sep = "")
		dir.create(outpath, showWarnings = F)
		
		if(!file.exists(paste(outpath, "/TICtable.csv", sep = ""))){
			data_file = dir(inpath, full.names = T)
			data_name = dir(inpath)
			data_name = gsub("(.*?).mzXML","\\1",data_name)
			totalRT = NULL
			TICtable = NULL
			pb <- tkProgressBar("TIC Clustering-Please wait TIC data extracting", "0% done", 0, 100, 0, width = 500)
			cat("TIC Clustering-Please wait TIC data extracting \n0 ")
			for(x in 1:length(data_file)){
				a = openMSfile(data_file[x])
				
				hd = header(a)
				
				totalRT = cbind(totalRT, hd[,"retentionTime"])
				TICtable = cbind(TICtable, hd[,"totIonCurrent"])
				Sys.sleep(0.1)
				info <- sprintf("%d%%", round(100*x/length(data_file)))
				setTkProgressBar(pb, value = round(100*x/length(data_file)), sprintf("TIC Clustering-Please wait TIC data extracting (%s done)", info), paste(data_name[x], info, sep = " "))
				if(round(100*x/length(data_file))<100){
					cat(round(100*x/length(data_file)), " ", sep = "")
				}else{
					cat(round(100*x/length(data_file)), " \n", sep = "")
				}
			}
			TICtable = cbind(1:nrow(TICtable), apply(totalRT, 1, mean), TICtable)
			colnames(TICtable) = c("Scan_index", "Ret_Time", gsub("\\d{4,8}_(.*?)", "\\1", data_name))
			write.table(TICtable, file = paste(outpath, "/TICtable.csv", sep = ""), col.names = T, row.names = F, sep = ",", quote = F)
			setTkProgressBar(pb, value = 100, "TIC Clustering-Finished TIC data extracting and start clustering.", "Finished 100%")
			Sys.sleep(1)
			close(pb)
		}else{
			TICtable = read.table(paste(outpath, "/TICtable.csv", sep = ""), sep = ",", header = T, quote = "\"", as.is = T, fill = T)
		}
		rt_start = as.numeric(tclvalue(RT_start))
		rt_end = as.numeric(tclvalue(RT_end))
		method = c("complete", "average", "ward")[as.numeric(c(tclvalue(method_c.val), tclvalue(method_a.val), tclvalue(method_w.val)))*c(1,2,3)]
		
		RT_idx = which(TICtable[,2]>rt_start & TICtable[,2]<rt_end)
		data_TIC = TICtable[RT_idx,-c(1,2)]
		ss= colnames(data_TIC)
		sample.info = data.frame(Sample = ss, Group = gsub("_[0-9]", "", ss), stringsAsFactors = F)

		r.cluster(data_TIC, r = rep(1, nrow(data_TIC)), cut.r = Inf, num.replicate = as.numeric(names(which.max(table(table(sample.info$Group))))), sample.info.temp = sample.info, c.path = outpath, index.plot = T, index.plot.sh = F, index.correlation = T, index.output = F, height_method = 1, method = method, QC = F)
		paste(toupper(substring(method, 1,1)), substring(method, 2), sep = "")
		file.rename(paste(outpath, toupper(substring(method, 1,1)), substring(method, 2), "/Clustering.png", sep = ""), paste(outpath, toupper(substring(method, 1,1)), substring(method, 2), "/Clustering_", rt_start, "-", rt_end, ".png", sep = ""))
		
		tkmessageBox(title="TIC Clustering",message="TIC Clustering is done.", icon = "info", type = "ok")
		cat("TIC Clustering - TIC Clustering is done.\n", sep = "")
		tkconfigure(tt,cursor="arrow")
	}
	onCancel <- function()
	{
		
		
		tkdestroy(dlg)
		tkfocus(tt)
	}
	fr <- tkframe(dlg)
	OK.but     <-tkbutton(fr,text=" Clustering ",command=onClust)
	Cancel.but <-tkbutton(fr,text="   Cancel   ",command=onCancel)
	tkgrid(tklabel(fr,text="               "), OK.but,tklabel(fr,text="                                      "), Cancel.but, tklabel(fr,text="               "))
	tkgrid(fr)
	tkgrid(tklabel(dlg, text = "", height = 0, font = fontIntro_para))

	tkfocus(dlg)
	tkbind(dlg, "<Destroy>", function() {tkgrab.release(dlg);tkfocus(tt)})
	
	
	tkwait.window(dlg)

	
}  


QC_interactive <- function(acc, r, r.range, metric = NULL, method, outpath, height_method, QCtable_temp){
	outpath <- paste(outpath, c("Pearson", "Spearman", "Kendall")[metric], sep="/")
	
	
	width_sh = width_acc = 480
	width_clust = 640
	width_list = 320
	height_list = height_clust = height_sh = height_acc = switch(ws, 360, 270, 240)
	fnt.lab = fnt.axis = c("Helvetica", switch(ws, "9", "7", "7"), "roman")
	fnt.method = c("Helvetica", switch(ws, "9", "7", "7"), "bold")
	fnt.probe = c("Helvetica", switch(ws, "8", "6", "6"), "roman")
	fnt.probeB = c("Helvetica", switch(ws, "8", "6", "6"), "bold")
	fnt.thresh = c("Helvetica", switch(ws, "9", "7", "7"), "bold", "italic")
	SQ <- sqrt(2)
	diam <- switch(ws, 2.9, 2.2, 2)
	
	
	default_method <- ifelse(any(method==2),2,method[1]) 
	default_r_choi <- which(!is.na(sapply(acc,function(x) x$accuracy[[default_method]])))[1]
	


	
	Point_acc <- function(x, y, pch, col, fill, diam) {
		switch(as.character(pch), `0` = Point_acc(x, y, 22, col, 
			fill = "", diam), `1` = Point_acc(x, y, 21, col, fill = "", 
			diam), `2` = Point_acc(x, y, 24, col, fill = "", diam), 
			`3` = {
				tkcreate(can, "line", x, y + SQ * diam, x, y - 
				  SQ * diam, fill = col)
				tkcreate(can, "line", x + SQ * diam, y, x - SQ * 
				  diam, y, fill = col)
			}, `4` = {
				tkcreate(can, "line", x - diam, y - diam, x + 
				  diam, y + diam, fill = col)
				tkcreate(can, "line", x - diam, y + diam, x + 
				  diam, y - diam, fill = col)
			}, `5` = Point_acc(x, y, 23, col, fill = "", diam), `6` = Point_acc(x, 
				y, 25, col, fill = "", diam), `7` = {
				Point_acc(x, y, 4, col, fill, diam)
				Point_acc(x, y, 0, col, fill, diam)
			}, `8` = {
				Point_acc(x, y, 3, col, fill, diam)
				Point_acc(x, y, 4, col, fill, diam)
			}, `9` = {
				Point_acc(x, y, 3, col, fill, diam)
				Point_acc(x, y, 5, col, fill, diam)
			}, `10` = {
				Point_acc(x, y, 3, col, fill, diam/SQ)
				Point_acc(x, y, 1, col, fill, diam)
			}, `11` = {
				Point_acc(x, y, 2, col, fill, diam)
				Point_acc(x, y, 6, col, fill, diam)
			}, `12` = {
				Point_acc(x, y, 3, col, fill, diam/SQ)
				Point_acc(x, y, 0, col, fill, diam)
			}, `13` = {
				Point_acc(x, y, 4, col, fill, diam)
				Point_acc(x, y, 1, col, fill, diam)
			}, `14` = {
				tkcreate(can, "line", x - diam, y - diam, x, 
				  y + diam, fill = col)
				tkcreate(can, "line", x + diam, y - diam, x, 
				  y + diam, fill = col)
				Point_acc(x, y, 0, col, fill, diam)
			}, `15` = Point_acc(x, y, 22, col = col, fill = col, 
				diam), `16` = Point_acc(x, y, 21, col = col, fill = col, 
				diam), `17` = Point_acc(x, y, 24, col = col, fill = col, 
				diam), `18` = Point_acc(x, y, 23, col = col, fill = col, 
				diam/SQ), `19` = Point_acc(x, y, 21, col = col, fill = col, 
				diam), `20` = Point_acc(x, y, 21, col = col, fill = col, 
				diam/2), `21` = tkcreate(can, "oval", x - diam-0.5, 
				y - diam-0.5, x + diam+0.5, y + diam+0.5, outline = col, 
				fill = fill, activefill = col), `22` = tkcreate(can, "rectangle", 
				x - diam-0.5, y - diam-0.5, x + diam+0.5, y + diam+0.5, outline = col, 
				fill = fill, activefill = col), `23` = tkcreate(can, "polygon", 
				x, y + SQ * diam, x + SQ * diam, y, x, y - SQ * 
				  diam, x - SQ * diam, y, outline = col, fill = fill, activefill = col), 
			`24` = tkcreate(can, "polygon", x, y - SQ * diam, 
				x + sqrt(6)/2 * diam, y + SQ/2 * diam, x - sqrt(6)/2 * 
				  diam, y + SQ/2 * diam, outline = col, fill = fill, activefill = col), 
			`25` = tkcreate(can, "polygon", x, y + SQ * diam, 
				x + sqrt(6)/2 * diam, y - SQ/2 * diam, x - sqrt(6)/2 * 
				  diam, y - SQ/2 * diam, outline = col, fill = fill, activefill = col), 
			o = Point_acc(x, y, 1, col, fill, diam),
			{tkcreate(can, "text", x, y, text = as.character(pch), 
				  fill = col)
				Point_acc(x, y, 21, col = "", fill = "", diam)
			})
    }
	
	Point_sh <- function(x, y, pch, col, fill, diam) {
		switch(as.character(pch), `0` = Point_sh(x, y, 22, col, 
			fill = "", diam), `1` = Point_sh(x, y, 21, col, fill = "", 
			diam), `2` = Point_sh(x, y, 24, col, fill = "", diam), 
			`3` = {
				tkcreate(can_sh, "line", x, y + SQ * diam, x, y - 
				  SQ * diam, fill = col)
				tkcreate(can_sh, "line", x + SQ * diam, y, x - SQ * 
				  diam, y, fill = col)
			}, `4` = {
				tkcreate(can_sh, "line", x - diam, y - diam, x + 
				  diam, y + diam, fill = col)
				tkcreate(can_sh, "line", x - diam, y + diam, x + 
				  diam, y - diam, fill = col)
			}, `5` = Point_sh(x, y, 23, col, fill = "", diam), `6` = Point_sh(x, 
				y, 25, col, fill = "", diam), `7` = {
				Point_sh(x, y, 4, col, fill, diam)
				Point_sh(x, y, 0, col, fill, diam)
			}, `8` = {
				Point_sh(x, y, 3, col, fill, diam)
				Point_sh(x, y, 4, col, fill, diam)
			}, `9` = {
				Point_sh(x, y, 3, col, fill, diam)
				Point_sh(x, y, 5, col, fill, diam)
			}, `10` = {
				Point_sh(x, y, 3, col, fill, diam/SQ)
				Point_sh(x, y, 1, col, fill, diam)
			}, `11` = {
				Point_sh(x, y, 2, col, fill, diam)
				Point_sh(x, y, 6, col, fill, diam)
			}, `12` = {
				Point_sh(x, y, 3, col, fill, diam/SQ)
				Point_sh(x, y, 0, col, fill, diam)
			}, `13` = {
				Point_sh(x, y, 4, col, fill, diam)
				Point_sh(x, y, 1, col, fill, diam)
			}, `14` = {
				tkcreate(can_sh, "line", x - diam, y - diam, x, 
				  y + diam, fill = col)
				tkcreate(can_sh, "line", x + diam, y - diam, x, 
				  y + diam, fill = col)
				Point_sh(x, y, 0, col, fill, diam)
			}, `15` = Point_sh(x, y, 22, col = col, fill = col, 
				diam), `16` = Point_sh(x, y, 21, col = col, fill = col, 
				diam), `17` = Point_sh(x, y, 24, col = col, fill = col, 
				diam), `18` = Point_sh(x, y, 23, col = col, fill = col, 
				diam/SQ), `19` = Point_sh(x, y, 21, col = col, fill = col, 
				diam), `20` = Point_sh(x, y, 21, col = col, fill = col, 
				diam/2), `21` = tkcreate(can_sh, "oval", x - diam-0.5, 
				y - diam-0.5, x + diam+0.5, y + diam+0.5, outline = col, 
				fill = fill), `22` = tkcreate(can_sh, "rectangle", 
				x - diam-0.5, y - diam-0.5, x + diam+0.5, y + diam+0.5, outline = col, 
				fill = fill), `23` = tkcreate(can_sh, "polygon", 
				x, y + SQ * diam, x + SQ * diam, y, x, y - SQ * 
				  diam, x - SQ * diam, y, outline = col, fill = fill), 
			`24` = tkcreate(can_sh, "polygon", x, y - SQ * diam, 
				x + sqrt(6)/2 * diam, y + SQ/2 * diam, x - sqrt(6)/2 * 
				  diam, y + SQ/2 * diam, outline = col, fill = fill), 
			`25` = tkcreate(can_sh, "polygon", x, y + SQ * diam, 
				x + sqrt(6)/2 * diam, y - SQ/2 * diam, x - sqrt(6)/2 * 
				  diam, y - SQ/2 * diam, outline = col, fill = fill), 
			o = Point_sh(x, y, 1, col, fill, diam),
			{tkcreate(can_sh, "text", x, y, text = as.character(pch), 
				  fill = col)
				Point_sh(x, y, 21, col = "", fill = "", diam)
			})
    }
	

	
	
	usr2xy <- function(row, xincr, yincr, mar, xy0) {
		x <- (row[1] - xy0[1]) * xincr + mar[2]
		y <- (xy0[2] - row[2]) * yincr + mar[3]
		c(x, y)
	}
	
	xy2usr <- function(item, xrange, yrange, xincr, yincr, mar) {
		xy <- as.numeric(tkcoords(can, item))
		x <- xy[1]
		y <- xy[2]
		x <- xrange[1] + (x - mar[2])/xincr
		y <- yrange[2] - (y - mar[3])/yincr
		c(x, y)
	}
	x2usr <- function(xcan, xrange, xincr, mar) {
		xrange[1] + (xcan - mar[2])/xincr
	}
	y2usr <- function(ycan, yrange, yincr, mar) {
		yrange[2] - (ycan - mar[3])/yincr
	}
	
	qc <- tktoplevel(); if(isIcon) tk2ico.set(qc,icon)
	tktitle(qc) <- "Quality Control - Output viewer"
	
	qc_txt <- tkframe(qc)
	tkpack(qc_txt, side = "top", fill = "x")
	lab_method <- tklabel(qc_txt, text = paste("   Clustering method : ", sep = ""))
	tkgrid(lab_method, sticky = "w")
	lab_r_thresh <- tklabel(qc_txt, text = paste("   r-value cut-off : ", sep = ""))
	tkgrid(lab_r_thresh, sticky = "w")
	lab_height <- tklabel(qc_txt, text = paste("   s-value threshold () : ", sep = ""))
	tkgrid(lab_height, sticky = "w")

	
	qc_up <- tkframe(qc)
	tkpack(qc_up, side = "top")
	fr_acc <- tkframe(qc_up, width = width_acc, height = height_acc) 
	accuracy = matrix(sapply(acc, function(x) unlist(x$accuracy)), nrow = length(r.range), byrow = T)
	accuracy = cbind(1:length(r.range), accuracy)
	rownames(accuracy) = as.character(r.range)
	
	
	
	pcol_acc = lcol_acc = c("red", "blue", "green")[method]
	tcol_acc = rep("black", length(method))
	tcol_probe = c("gray47", "gray87")[(1:length(r.range))%%2 + 1]
	pbg_acc = c("red", "blue", "green")
	pch_acc = c(21, 22, 24)[method]
	xlim <- range(accuracy[, 1], na.rm = TRUE)
	ylim <- range(accuracy[,-1], na.rm = TRUE)*c(0.8, 1)
	xpretty <- pretty(xlim)
	ypretty <- pretty(ylim)
	xrange <- c(-0.04, 0.04) * diff(xlim) + xlim
	xpretty <- xpretty[xpretty >= xrange[1] & xpretty <= xrange[2]]
	yrange <- c(-0.04, 0.04) * diff(ylim) + ylim
	ypretty <- ypretty[ypretty >= yrange[1] & ypretty <= yrange[2]]
	rpix <- 19.2
	mar <- round(c(2, 3.3, 1.5, 0.5) * rpix)
	xusr <- width_acc - mar[2] - mar[4]
	xincr <- xusr/diff(xrange)
	yusr <- height_acc -mar[1] - mar[3]
	yincr <- yusr/diff(yrange)
	xy0 <- c(xrange[1], yrange[2])
	
	
	can <- tkcanvas(fr_acc, relief = "sunken", width = width_acc, height = height_acc, scrollregion = c(0, 0, width_acc, height_acc))
	tkpack(can, side = "left", fill = "x")
	x0 <- usr2xy(c(xrange[1], yrange[1]), xincr, yincr, mar, xy0)
	x1 <- usr2xy(c(xrange[2], yrange[2]), xincr, yincr, mar, xy0)
	tkcreate(can, "rectangle", x0[1], x0[2], x1[1], x1[2], outline = "black", width = 1, fill = "white")
	tl <- rpix*0.25
	
	
	
	
	
	tmp <- rownames(accuracy)
	for (i in 1:length(tmp)) {  
		xx <- usr2xy(c(accuracy[i, 1], yrange[1]), xincr, yincr, mar, xy0)
		tkcreate(can, "line", xx[1], xx[2], xx[1], xx[2] + tl, fill = "black")
		tkcreate(can, "text", xx[1], xx[2] + rpix * 0.5, 
			anchor = "n", text = as.character(tmp[i]), fill = "black", 
			font = fnt.axis)
		item <- tkcreate(can, "text", xx[1], xx[2] - rpix * 0.4, 
			text = as.character(acc[[i]]$num.probe), fill = tcol_probe[i], 
			font = fnt.probe)
		tkaddtag(can, paste("text-", i, sep = ""), "withtag", item)
		
			
			
	}
	
	xx <- usr2xy(c(mean(xrange), yrange[1]), xincr, yincr, mar, xy0)
	tkcreate(can, "text", xx[1], xx[2] + rpix * 1.13, text = "r-value cut-off", 
		fill = "black", anchor = "n", font = fnt.lab)
	
	tmp <- ypretty
	for (i in 1:length(tmp)) {    
		
		
		
		yy <- usr2xy(c(xrange[1], tmp[i]), xincr, yincr, mar, xy0)
		tkcreate(can, "line", yy[1], yy[2], yy[1] - tl, yy[2], fill = "black")
		tkcreate(can, "text", yy[1] - rpix * 0.5, yy[2], 
			anchor = "e", text = as.character(tmp[i]), fill = "black", 
			font = fnt.axis)
	} 

	yy <- usr2xy(c(xrange[1], yrange[2]), xincr, yincr, mar, xy0)
	tkcreate(can, "text", yy[1] - rpix * 0.3, yy[2], text = "Accuracy", fill = "black", anchor = "se", font = fnt.lab)
	
	
	
	for (j in 1:length(method)) { 
		for (i in 1:(nrow(accuracy)-1)) {
			if(any(is.na(accuracy[i, c(1, j+1)]))) next
			xy <- usr2xy(accuracy[i, c(1, j+1)], xincr, yincr, mar, xy0)
			
				
			if(any(is.na(accuracy[i+1, c(1, j+1)]))) next
			tmp <- usr2xy(accuracy[i+1, c(1, j+1)], xincr, yincr, mar, xy0)
			item <- tkcreate(can, "line", xy[1], xy[2], tmp[1], tmp[2], fill = pcol_acc[j])
			
			
			
				
			tkaddtag(can, paste("line-", c("complete", "average", "ward")[method[j]], sep = ""), "withtag", item)
			
			
			
			
		}
		
		xy <- usr2xy(c(xrange[1], yrange[1]), xincr, yincr, mar, xy0)
		tkcreate(can, "line", xy[1] + tl, xy[2] - rpix * c(1.2, 2, 2.8)[j], xy[1] + 4 * tl, xy[2] - rpix * c(1.2, 2, 2.8)[j], fill = pcol_acc[j])
		item <- tkcreate(can, "text", xy[1] + 5*tl + rpix * 0.5, xy[2] - rpix * c(1.2, 2, 2.8)[j], 
			anchor = "w", text = c("Complete linkage", "Average linkage", "Ward's method")[method[j]], fill = pcol_acc[j], 
			font = fnt.method)
		tkaddtag(can, paste("method-", c("complete", "average", "ward")[method[j]], sep = ""), "withtag", item)
	}
	for (j in 1:length(method)) { 
		for (i in 1:nrow(accuracy)) {
			if(any(is.na(accuracy[i, c(1, j+1)]))) next

			xy <- usr2xy(accuracy[i, c(1, j+1)], xincr, yincr, mar, xy0)
			if(j==1 & i==1){
				item <- Point_acc(xy[1], xy[2], pch = pch_acc[j], col = pcol_acc[j], 
					fill = pcol_acc[j], diam = diam)
			}else{
				item <- Point_acc(xy[1], xy[2], pch = pch_acc[j], col = pcol_acc[j], 
					fill = "white", diam = diam)
			}
			
			
				
				
			
			
			
			
				
			tkaddtag(can, paste("point-", c("complete", "average", "ward")[method[j]], sep = ""), "withtag", item)
			tkaddtag(can, paste("clust-", c("complete", "average", "ward")[method[j]], sep = ""), "withtag", item)
			tkaddtag(can, paste("list-", c("complete", "average", "ward")[method[j]], sep = ""), "withtag", item)
			
			
			
			
			
		}
	}
	pEnter_acc <- function() {
		tkconfigure(qc,cursor="hand2")
		
		
			
			
			
				
		
			
				
		
		
		
		
		
		default_r_choi <<- round(x2usr(mean(as.numeric(tkbbox(can, "current"))[c(1,3)]), xrange, xincr, mar))
		default_method <<- which(c("point-complete", "point-average", "point-ward")==as.character(tkgettags(can, "current"))[1])
		tkitemconfigure(can, paste("text-", default_r_choi, sep = ""), fill = "magenta", font = fnt.probeB)
	}
	pLeave_acc <- function() {
		tkconfigure(qc,cursor="arrow")
		
		
		tkitemconfigure(can, paste("text-", default_r_choi, sep = ""), fill = tcol_probe[default_r_choi], font = fnt.probe)
	}
	
	
	
	fr_sh <- tkframe(qc_up, width = width_sh, height = height_sh)

	
	tkgrid(fr_acc, fr_sh)
	
	
	sh_interactive <- function(){
		if(exists("can_sh")) tkpack.forget(can_sh)
		
		
		tkitemconfigure(can, "point-complete", fill = "white")
		tkitemconfigure(can, "point-average", fill = "white")
		tkitemconfigure(can, "point-ward", fill = "white")
		
		tkitemconfigure(can, "current", fill = pbg_acc[default_method])
		
		tkconfigure(lab_r_thresh, text = paste("   r-value cut-off : ", r.range[default_r_choi], sep = ""))
		tkconfigure(lab_method, text = paste("   Clustering method : ", c("Complete linkage", "Average linkage", "Ward's method")[default_method], sep = ""))
		tkconfigure(lab_height, text = paste("   s-value threshold (", c("Parametric", "Nonparametric", "Bootstrap")[height_method], ") : ", round(acc[[default_r_choi]]$height_thresh[[c("complete", "average", "ward")[default_method]]], 3), sep = ""))
		
		
		
		
		step_height = cbind(acc[[default_r_choi]]$steps[[c("complete", "average", "ward")[default_method]]], acc[[default_r_choi]]$height[[c("complete", "average", "ward")[default_method]]])
		height_thresh = acc[[default_r_choi]]$height_thresh[[c("complete", "average", "ward")[default_method]]]
		index_bad = NA
		index_h2 = NA
		if(length(acc[[default_r_choi]]$index.bad[[c("complete", "average", "ward")[default_method]]]))
		{
			index_bad = acc[[default_r_choi]]$index.bad[[c("complete", "average", "ward")[default_method]]]
			labels1 = acc[[default_r_choi]]$sb.bad[[c("complete", "average", "ward")[default_method]]]
			labels1 = labels1[order(labels1)]
		}
		if(length(acc[[default_r_choi]]$index.h2[[c("complete", "average", "ward")[default_method]]]))
		{
			index_h2 = acc[[default_r_choi]]$index.h2[[c("complete", "average", "ward")[default_method]]]
			labels2 = acc[[default_r_choi]]$sb.bad2[[c("complete", "average", "ward")[default_method]]]
			labels2 = labels2[order(labels2)]
		}

		
		
		
		
		
		
		
		
		
		
		pcol_sh = "black"
		tcol_sh = "gray38"
		pbg_sh = acc[[default_r_choi]]$backcol[[c("complete", "average", "ward")[default_method]]]
		pbg_sh[pbg_sh==1] = "black"
		pbg_sh[pbg_sh==2] = "red"
		pbg_sh[pbg_sh==3] = "green"
		
		pch_sh = acc[[default_r_choi]]$point[[c("complete", "average", "ward")[default_method]]]
		xlim <- range(step_height[, 1], na.rm = TRUE)
		if(diff(xlim)<2) xlim <- c(xlim[1]-1, xlim[2]+1)
		ylim <- c(range(step_height[, 2], na.rm = TRUE)[1], max(range(step_height[, 2], na.rm = TRUE)[2], height_thresh))
		xpretty <- pretty(xlim)
		ypretty <- pretty(ylim)
		xrange <- c(-0.04, 0.04) * diff(xlim) + xlim
		xpretty <- xpretty[xpretty >= xrange[1] & xpretty <= xrange[2]]
		yrange <- c(-0.04, 0.04) * diff(ylim) + ylim
		ypretty <- ypretty[ypretty >= yrange[1] & ypretty <= yrange[2]]
		rpix <- 19.2
		mar <- round(c(2, 2.5, 1.5, 0.5) * rpix)
		xusr <- width_sh - mar[2] - mar[4]
		xincr <- xusr/diff(xrange)
		yusr <- height_sh -mar[1] - mar[3]
		yincr <- yusr/diff(yrange)
		xy0 <- c(xrange[1], yrange[2])
		
		can_sh <<- tkcanvas(fr_sh, relief = "sunken", width = width_sh, height = height_sh, scrollregion = c(0, 0, width_sh, height_sh))
		tkpack(can_sh, side = "left")
		x0 <- usr2xy(c(xrange[1], yrange[1]), xincr, yincr, mar, xy0)
		x1 <- usr2xy(c(xrange[2], yrange[2]), xincr, yincr, mar, xy0)
		tkcreate(can_sh, "rectangle", x0[1], x0[2], x1[1], x1[2], outline = "black", width = 1, fill = "white")
		x0_thresh <- usr2xy(c(xrange[1], height_thresh), xincr, yincr, mar, xy0)
		x1_thresh <- usr2xy(c(xrange[2], height_thresh), xincr, yincr, mar, xy0)
		tkcreate(can_sh, "line", x0_thresh[1], x0_thresh[2], x1_thresh[1], x1_thresh[2], width = 1, fill = "green")
		tkcreate(can_sh, "text", x0_thresh[1], x0_thresh[2], anchor = "sw", text = paste("s-value=", as.character(round(height_thresh, 3)), sep = ""), fill = "green", font = fnt.thresh)
		tl <- rpix*0.25
		tmp <- xpretty
		for (i in 1:length(tmp)) {  
			xx <- usr2xy(c(tmp[i], yrange[1]), xincr, yincr, mar, xy0)
			tkcreate(can_sh, "line", xx[1], xx[2], xx[1], xx[2] + tl, fill = "black")
			tkcreate(can_sh, "text", xx[1], xx[2] + rpix * 0.5, 
				anchor = "n", text = as.character(tmp[i]), fill = "black", 
				font = fnt.axis)
		}
		
		xx <- usr2xy(c(mean(xrange), yrange[1]), xincr, yincr, mar, xy0)
		tkcreate(can_sh, "text", xx[1], xx[2] + rpix * 1, text = "Steps", 
			fill = "black", anchor = "n", font = fnt.lab)
		
		tmp <- ypretty
		for (i in 1:length(tmp)) {    
			
			
			
			yy <- usr2xy(c(xrange[1], tmp[i]), xincr, yincr, mar, xy0)
			tkcreate(can_sh, "line", yy[1], yy[2], yy[1] - tl, yy[2], fill = "black")
			tkcreate(can_sh, "text", yy[1] - rpix * 0.5, yy[2], 
				anchor = "e", text = as.character(tmp[i]), fill = "black", 
				font = fnt.axis)
		} 

		yy <- usr2xy(c(xrange[1], yrange[2]), xincr, yincr, mar, xy0)
		tkcreate(can_sh, "text", yy[1] - rpix * 0.3, yy[2], text = "h", fill = "black", anchor = "se", font = fnt.lab)
		
		temp_pch = c(21, 24, 21, 21, 21)
		temp_col = c("black", "black", "", "", "")
		temp_pbg = c("white", "white", "black", "green", "red")
		temp_txt = c(paste("h<", round(height_thresh,3), sep = ""), paste("h>", round(height_thresh,3), sep = ""), "perfectly-clustered", "well-clustered", "poorly-clustered")
		for(i in 1:5){
			xy <- usr2xy(c(xrange[1], yrange[2]-diff(yrange)*0.045*i), xincr, yincr, mar, xy0)
			Point_sh(xy[1] + rpix*0.5, xy[2], pch = temp_pch[i], col = temp_col[i], fill = temp_pbg[i], diam = diam)
			tkcreate(can_sh, "text", xy[1] + rpix , xy[2], 
				anchor = "w", text = temp_txt[i], fill = "black", 
				font = fnt.method)
		}
		
		pola_sh <- tclArray()
		labtext_sh <- tclArray()
		id_sh <- tclArray()
		fnt <- c("Helvetica", "9", "roman")
		
		for (i in 1:nrow(step_height)) {
			if(i %in% c(index_bad, index_h2)) next
			xy <- usr2xy(step_height[i,], xincr, yincr, mar, xy0)
			item <- Point_sh(xy[1], xy[2], pch = pch_sh[i], col = pcol_sh, 
				fill = pbg_sh[i], diam = diam)
			
			tkaddtag(can_sh, "point", "withtag", item)
		}
		
		if(all(!is.na(index_bad))){
			temp = t(cbind(matrix(step_height[index_bad,], nrow = length(index_bad)), pch_sh[index_bad], pbg_sh[index_bad]))
			temp = unique(t(subset(temp, select = order(temp[1,], temp[2,]))))
			temp = matrix(temp[order(rownames(temp)),], ncol = 4)
			for(i in 1:nrow(temp)){
				xy <- usr2xy(as.numeric(temp[i,1:2]), xincr, yincr, mar, xy0)
				item <- Point_sh(xy[1], xy[2], pch = as.numeric(temp[i,3]), col = pcol_sh, 
					fill = temp[i,4], diam = diam)
				
				tkaddtag(can_sh, "point", "withtag", item)
				lab <- tkcreate(can_sh, "text", xy[1], xy[2]-rpix*0.5, 
					text = labels1[i], fill = tcol_sh, font = fnt)
				tkaddtag(can_sh, "label", "withtag", lab)
				pola_sh[[lab]] <- item
				labtext_sh[[lab]] <- labels1[i]
				id_sh[[lab]] <- i
			}
		}
		if(all(!is.na(index_h2))){
			temp = t(cbind(matrix(step_height[index_h2,], nrow = length(index_h2)), pch_sh[index_h2], pbg_sh[index_h2]))
			temp = unique(t(subset(temp, select = order(temp[1,], temp[2,]))))
			temp = matrix(temp[order(rownames(temp)),], ncol = 4)
			for(i in 1:nrow(temp)){
				xy <- usr2xy(as.numeric(temp[i,1:2]), xincr, yincr, mar, xy0)
				item <- Point_sh(xy[1], xy[2], pch = as.numeric(temp[i,3]), col = pcol_sh, 
					fill = temp[i,4], diam = diam)
				
				tkaddtag(can_sh, "point", "withtag", item)
				lab <- tkcreate(can_sh, "text", xy[1], xy[2]-rpix*0.5, 
					text = labels2[i], fill = tcol_sh, font = fnt)
				tkaddtag(can_sh, "label", "withtag", lab)
				pola_sh[[lab]] <- item
				labtext_sh[[lab]] <- labels2[i]
				id_sh[[lab]] <- i
			}
		}
		
		pEnter_sh <- function() {
			tkdelete(can_sh, "box")
			hbox <- tkcreate(can_sh, "rectangle", tkbbox(can_sh, "current"), 
				outline = "red", fill = "yellow")
			tkaddtag(can_sh, "box", "withtag", hbox)
			tkitemraise(can_sh, "current")
		}
		pLeave_sh <- function() {
			tkdelete(can_sh, "box")
		}
		pDown_sh <- function(x, y) {
			x <- as.numeric(x)
			y <- as.numeric(y)
			tkdtag(can_sh, "selected")
			tkaddtag(can_sh, "selected", "withtag", "current")
			tkitemraise(can_sh, "current")
			p <- as.numeric(tkbbox(can_sh, pola_sh[[tkfind(can_sh, "withtag", 
				"current")]]))
			.pX <<- (p[1] + p[3])/2
			.pY <<- (p[2] + p[4])/2
			.lastX <<- x
			.lastY <<- y
		}
		pMove_sh <- function(x, y) {
			
			x <- as.numeric(x)
			if(x<mar[2]) x <- mar[2]
			if(x>(width_sh - mar[4])) x <- width_sh - mar[4]
			y <- as.numeric(y)
			if(y<mar[3]) y <- mar[3]
			if(y>(height_sh - mar[1])) y <- height_sh - mar[1]
			tkmove(can_sh, "selected", x - .lastX, y - .lastY)
			tkdelete(can_sh, "ptr")
			tkdelete(can_sh, "box")
			.lastX <<- x
			.lastY <<- y
			xadj <- as.numeric(tkcanvasx(can_sh, 0))
			yadj <- as.numeric(tkcanvasy(can_sh, 0))
			hbox <- tkcreate(can_sh, "rectangle", tkbbox(can_sh, "selected"), 
				outline = "red")
			tkaddtag(can_sh, "box", "withtag", hbox)
			conn <- tkcreate(can_sh, "line", .lastX + xadj, .lastY + 
				yadj, .pX, .pY, fill = "red")
			tkaddtag(can_sh, "ptr", "withtag", conn)
		}
		.lastX <- 0
		.lastY <- 0
		.pX <- 0
		.pY <- 0
		tkitembind(can_sh, "label", "<Any-Enter>", pEnter_sh)
		tkitembind(can_sh, "label", "<Any-Leave>", pLeave_sh)
		tkitembind(can_sh, "label", "<1>", pDown_sh)
		tkitembind(can_sh, "label", "<ButtonRelease-1>", function() {
			tkdtag(can_sh, "selected")
			tkdelete(can_sh, "ptr")
		})
		tkitembind(can_sh, "label", "<B1-Motion>", pMove_sh)
	}
	sh_interactive()

	qc_down <- tkframe(qc)
	tkpack(qc_down, side = "top")
	
	fr_clust <- tkframe(qc_down, width = width_clust, height = height_clust)
	fr_list <- tkframe(qc_down, width = width_list, height = height_list)
	
	tkgrid(fr_clust, tklabel(qc_down, text = "", width = 2), fr_list)
	
	clust_interactive <- function(){
		if(exists("can_clust")) tkpack.forget(can_clust)
		
		
		img_clust <- tclVar()

		tkimage.create("photo", img_clust, file = paste(outpath, c("/Complete", "/Average", "/Ward")[default_method], "/Clustering_r_", r.range[default_r_choi], ".png", sep = ""))
		
		img2_clust <- tclVar()
		tkimage.create("photo", img2_clust)
		tcl(img2_clust, "copy", img_clust, zoom = 2)

		img_tmp <- readPNG(paste(outpath, c("/Complete", "/Average", "/Ward")[default_method], "/Clustering_r_", r.range[default_r_choi], ".png", sep = ""), T)
		png(paste(Sys.getenv("subfunc"), "/zoom2.png", sep = ""), width = 1280, height = 720)
		grid.raster(img_tmp, interpolate = T)
		dev.off()
		
		png(paste(Sys.getenv("subfunc"), "/zoom1.png", sep = ""), width = 640, height = 360)
		grid.raster(img_tmp, interpolate = T)
		dev.off()
		
		
		img3_clust <- tclVar()
		tkimage.create("photo", img3_clust, file = paste(Sys.getenv("subfunc"), "/zoom2.png", sep = ""))
		
		
		
		img4_clust <- tclVar()
		tkimage.create("photo", img4_clust, file = paste(Sys.getenv("subfunc"), "/zoom1.png", sep = ""))
		
		
		can_clust <<- tkcanvas(fr_clust, width = width_clust, height = height_clust, scrollregion = c(0, 0, width_clust, height_clust))
		tkpack(can_clust, side = "left")
		i1 <- tclVar()
		tkimage.create("photo", i1)
		tcl(i1, "copy", img4_clust, "-from", 0,0, width_clust, height_clust)

		img_item <- tkcreate(can_clust, "image", width_clust/2, height_clust/2, image = i1)
		
		pDown_clust <- function(x, y) {
			x <- as.numeric(x)
			y <- as.numeric(y)
			.lastx <<- x
			.lasty <<- y
		}
		
		width_tmp <<- 640
		height_tmp <<- 360
		zoom <<- 1
		
		pMove_clust <- function(x, y) {
			x <- as.numeric(x)
			y <- as.numeric(y)
			
			img_tmp = switch(zoom, img4_clust, img3_clust, img_clust, img2_clust)
			tcl(i1, "copy", img_tmp, "-from", min(max(.lastX+.lastx-x, 0), width_tmp-width_clust), min(max(.lastY+.lasty-y, 0), height_tmp-height_clust), min(.lastX+.lastx-x + width_clust, width_tmp), min(.lastY+.lasty-y + height_clust, height_tmp))
			
			.lastX <<- min(max(.lastX+.lastx-x, 0), width_tmp-width_clust)
			.lastY <<- min(max(.lastY+.lasty-y, 0), height_tmp-height_clust)
			.lastx <<- x
			.lasty <<- y
		}
		
			
			
			
			
			
			
			
			
			
		
		pZoom_clust <- function(){
			if(zoom < 4){
				zoom <<- zoom + 1
				width_tmp <<- switch(zoom-1, 1280, 1920, 3840)
				height_tmp <<- switch(zoom-1, 720, 1080, 2160)
				tmp = switch(zoom-1, 2, 1.5, 2)
				img_tmp = switch(zoom-1, img3_clust, img_clust, img2_clust)
				tcl(i1, "copy", img_tmp, "-from", round(tmp*.lastX), round(tmp*.lastY), min(width_tmp, round(tmp*.lastX)+width_clust), min(height_tmp, round(tmp*.lastY)+height_clust))
				tkitembind(can_clust, img_item, "<B1-Motion>", pMove_clust)
			}
		}
		pOrin_clust <- function(){
			if(zoom > 1){
				zoom <<- zoom - 1
				width_tmp <<- switch(zoom, 640, 1280, 1920)
				height_tmp <<- switch(zoom, 360, 720, 1080)
				tmp = switch(zoom, 0.5, 2/3, 0.5)
				img_tmp = switch(zoom, img4_clust, img3_clust, img_clust)
				tcl(i1, "copy", img_tmp, "-from", min(.lastX, width_tmp-width_clust), min(.lastY, height_tmp-height_clust), min(width_tmp, .lastX+width_clust), min(height_tmp, .lastY+height_clust))
				tkitembind(can_clust, img_item, "<B1-Motion>", pMove_clust)
			}
		}
		.lastX <<- 0
		.lastY <<- 0
		tkitembind(can_clust, img_item, "<1>", pDown_clust)
		tkitembind(can_clust, img_item, "<B1-Motion>", pMove_clust)
		tkitembind(can_clust, img_item, "<Button-3>", pZoom_clust)
		tkitembind(can_clust, img_item, "<Double-Button-1>", pOrin_clust)
	}
	clust_interactive()
	
	
	list_interactive = function(){
		if(exists("fr_list_y")){
			tkpack.forget(fr_list_y)
			tkpack.forget(fr_list_x)
		}
		QC_list = read.table(paste(outpath, c("/Complete", "/Average", "/Ward")[default_method], "/QC_r_", r.range[default_r_choi], ".csv", sep = ""), skip = 3, header = F, as.is = T, sep = ",")
		QC_list = cbind(c(0, 1:(nrow(QC_list)-1)), QC_list)
		QC_list[1,] = c("Idx", "SampleID", "Remove", "Quality", "h")
		list_array <- tclArray()
		for(i in 1:nrow(QC_list))
			for(j in 1:ncol(QC_list))
				list_array[[i-1,j-1]] <- QC_list[i,j]
		
		fr_list_y <<- tkframe(fr_list)
		tkpack(fr_list_y)
		yscr_list <<- tkscrollbar(fr_list_y, repeatinterval = 20, orient="vertical", command = function(...) tkyview(table_list,...))
		fr_list_x <<- tkframe(fr_list)
		tkpack(fr_list_x, fill = "x")
		xscr_list <<- tkscrollbar(fr_list_x, repeatinterval = 1, orient="horizontal", command = function(...) tkxview(table_list,...))
		table_list <<- tkwidget(fr_list_y, "table", variable = list_array, rows=nrow(QC_list), cols=5, titlerows = 1, background = "white", yscrollcommand = function(...) tkset(yscr_list, ...), xscrollcommand = function(...) tkset(xscr_list, ...), resizeborders="col", bordercursor = "sb_h_double_arrow", maxwidth = width_list-10, height = switch(ws, 20, 15, 13), padx = 10)
		
		
		
		tkpack(table_list, side = "left")
		tkpack(yscr_list, side = "left", fill = "y")
		tkpack(xscr_list, fill = "x")
		
		poor_idx = which(QC_list[,4]=="poor")
		well_idx = which(QC_list[,4]=="well")
		eliminate_idx = which(QC_list[,3]=="Yes")
		tcl(.Tk.ID(table_list), "tag", "configure", "poor", bg = "red", fg = "")
		tcl(.Tk.ID(table_list), "tag", "configure", "well", bg = "green", fg = "")
		tcl(.Tk.ID(table_list), "tag", "configure", "eliminate", bg = "yellow", fg = "")
		if(length(poor_idx)){
			for(i in 1:length(poor_idx)){
				tcl(.Tk.ID(table_list),"tag","celltag","poor",paste(poor_idx[i]-1, "3", sep = ","))
			}
		}
		if(length(well_idx)){
			for(i in 1:length(well_idx)){
				tcl(.Tk.ID(table_list),"tag","celltag","well",paste(well_idx[i]-1, "3", sep = ","))
			}
		}
		if(length(eliminate_idx)){
			for(i in 1:length(eliminate_idx)){
				tcl(.Tk.ID(table_list),"tag","celltag","eliminate",paste(eliminate_idx[i]-1, "2", sep = ","))
			}
		}
		tkconfigure(table_list,selectmode="extended",rowseparator="\"\n\"",colseparator="\"\t\"", state = "disable")
	}

	list_interactive()

	
	tkpack(tklabel(qc, text = "   ", font = fontIntro_para, height = 0), side = "top")
	realign <- tkframe(qc)
	realign.but <- tkbutton(realign, text = "   Re-alignment   ", command = function(...){
		
		cut.r = r.range[default_r_choi]
		
			
			
		
			
		
		
		
		
		
		
		tkdestroy(qc)
		
		Align(ALtitle = "Peak Re-alignment and Annotation", afterQC = T, QCpath = paste(outpath, c("/Complete", "/Average", "/Ward")[default_method], "/QC_r_", cut.r, ".csv", sep = ""))
	})
	close.but <- tkbutton(realign, text = "   Close   ", command = function(){
		ReturnVal <- tkmessageBox(title = "Quality Control", message = "Do you really want to close without doing re-alignment?", icon = "question", type = "yesno", default = "yes")
		cat("QC Viewer - Do you really want to close without doing re-alignment?", tclvalue(ReturnVal), "\n")
		if(tclvalue(ReturnVal)=="yes"){
			cut.r = r.range[default_r_choi]
			remain_sample = read.table(paste(outpath, c("/Complete", "/Average", "/Ward")[default_method], "/QC_r_", cut.r, ".csv", sep = ""), sep = ",", header = T, quote = "\"", as.is = T, skip = 3)
			QCtable_temp = QCtable_temp[which(r<cut.r),c(1:3, sort(match(remain_sample[remain_sample[,2]=="No",1], colnames(QCtable_temp))))]
			write.table(QCtable_temp, file = paste(outpath, "/peakTable_postQC_S", ncol(QCtable_temp)-3, "_P", nrow(QCtable_temp), ".csv", sep = ""), col.names = T, row.names = F, quote = F, sep = ",")
			textAbuninput <<- tclVar(paste(outpath, "/peakTable_postQC_S", ncol(QCtable_temp)-3, "_P", nrow(QCtable_temp), ".csv", sep = ""))
			tkentryconfigure(peakMenu, 1, state = "active")
			tkentryconfigure(batchsubMenu, 0, state = "active")
			tkentryconfigure(batchsubMenu, 1, state = "active")
			tkentryconfigure(statMenu, 1, state = "active")
			tkentryconfigure(statMenu, 2, state = "active")
			tkdestroy(qc)
			tkmessageBox(title = "Quality Control", message = paste("Post-QC peak abundance data is saved as ''peakTable_postQC_S", ncol(QCtable_temp)-3, "_P", nrow(QCtable_temp), ".csv''.",sep=""), icon = "info", type = "ok") 
			cat(paste("QC Viewer - Post-QC peak abundance data is saved as peakTable_postQC_s", ncol(QCtable_temp)-3, "_P", nrow(QCtable_temp), ".csv.\n",sep="")) 
		}else if(tclvalue(ReturnVal)=="no"){
			tkfocus(qc)
		}
	})
	tkgrid(tklabel(realign, text = "      "), realign.but, tklabel(realign, text = "                                                   "), close.but, tklabel(realign, text = "      "))
	tkgrid(tklabel(realign, text = "   ", font = fontIntro_para, height = 0))
	tkpack(realign, side = "top")
	
	tkitembind(can, "point-complete", "<Any-Enter>", pEnter_acc)
	tkitembind(can, "point-complete", "<Any-Leave>", pLeave_acc)
	tkitembind(can, "point-average", "<Any-Enter>", pEnter_acc)
	tkitembind(can, "point-average", "<Any-Leave>", pLeave_acc)
	tkitembind(can, "point-ward", "<Any-Enter>", pEnter_acc)
	tkitembind(can, "point-ward", "<Any-Leave>", pLeave_acc)
	
	
	tkitembind(can, "point-complete", "<Button-1>", sh_interactive)
	tkitembind(can, "clust-complete", "<Button-1>", clust_interactive)
	tkitembind(can, "list-complete", "<Button-1>", list_interactive)
	tkitembind(can, "point-average", "<Button-1>", sh_interactive)
	tkitembind(can, "clust-average", "<Button-1>", clust_interactive)
	tkitembind(can, "list-average", "<Button-1>", list_interactive)
	tkitembind(can, "point-ward", "<Button-1>", sh_interactive)
	tkitembind(can, "clust-ward", "<Button-1>", clust_interactive)
	tkitembind(can, "list-ward", "<Button-1>", list_interactive)
	
	tkitembind(can, "method-complete", "<Button-1>", function(){
		tkitemconfigure(can, "point-complete", state = "hidden")
		tkitemconfigure(can, "line-complete", state = "hidden")
		tkitemconfigure(can, "method-complete", font = fnt.lab)
	})
	tkitembind(can, "method-complete", "<Double-Button-1>", function(...){
		tkitemconfigure(can, "point-complete", state = "normal")
		tkitemconfigure(can, "line-complete", state = "normal")
		tkitemconfigure(can, "method-complete", font = fnt.method)
	})
	tkitembind(can, "method-average", "<Button-1>", function(){
		tkitemconfigure(can, "point-average", state = "hidden")
		tkitemconfigure(can, "line-average", state = "hidden")
		tkitemconfigure(can, "method-average", font = fnt.lab)
	})
	tkitembind(can, "method-average", "<Double-Button-1>", function(...){
		tkitemconfigure(can, "point-average", state = "normal")
		tkitemconfigure(can, "line-average", state = "normal")
		tkitemconfigure(can, "method-average", font = fnt.method)
	})
	tkitembind(can, "method-ward", "<Button-1>", function(){
		tkitemconfigure(can, "point-ward", state = "hidden")
		tkitemconfigure(can, "line-ward", state = "hidden")
		tkitemconfigure(can, "method-ward", font = fnt.lab)
	})
	tkitembind(can, "method-ward", "<Double-Button-1>", function(...){
		tkitemconfigure(can, "point-ward", state = "normal")
		tkitemconfigure(can, "line-ward", state = "normal")
		tkitemconfigure(can, "method-ward", font = fnt.method)
	})
	
	
}  



QualityControl <- function()
{
	dlg <- tktoplevel(width = 800); if(isIcon) tk2ico.set(dlg,icon)
	tkwm.title(dlg, "Peak/Sample filtering")
	
	fr_input <- tkframe(dlg, width = 800)
	
	
	tkgrid(tklabel(fr_input, text = "", height = 0, font = fontIntro_para))
	tkgrid(tklabel(fr_input, text = paste("Input file: ", tclvalue(textAbuninput), "    ", sep = ""), padx = 11, wraplength = 600, justify = "left"), sticky = "w")
	
	
	
	
	tkgrid(tklabel(fr_input, text = paste("Output folder: ", tclvalue(textoutput), "    ", sep = ""), padx = 11, wraplength = 600, justify = "left"), sticky = "w")
	tkgrid(tklabel(fr_input, text = "", height = 0, font = fontIntro_para))
	tkgrid(fr_input, sticky = "w")
	
	
	
	
	
	
	
	
	
	
	
	fr_rtmz <- tkframe(dlg)
	rtmz_choice <- tclVar("1")
	rtmz_yes <- tkradiobutton(fr_rtmz, variable = rtmz_choice, value = "1", command = function(...){
		tkconfigure(rtmz1, state = "normal")
		if(tclvalue(rt.val)=="1"){
			tkconfigure(rtstart.entry, state = "normal")
			tkconfigure(rtend.entry, state = "normal")
		}else{
			tkconfigure(rtstart.entry, state = "disable")
			tkconfigure(rtend.entry, state = "disable")
		}
		tkconfigure(rtmz2, state = "normal")
		if(tclvalue(mz.val)=="1"){
			tkconfigure(mzstart.entry, state = "normal")
			tkconfigure(mzend.entry, state = "normal")
		}else{
			tkconfigure(mzstart.entry, state = "disable")
			tkconfigure(mzend.entry, state = "disable")
		}
	})
	rtmz_no <- tkradiobutton(fr_rtmz, variable = rtmz_choice, value = "2", command = function(...){
		tkconfigure(rtmz1, state = "disable")
		tkconfigure(rtstart.entry, state = "disable")
		tkconfigure(rtend.entry, state = "disable")
		tkconfigure(rtmz2, state = "disable")
		tkconfigure(mzstart.entry, state = "disable")
		tkconfigure(mzend.entry, state = "disable")
	})
	
	rt.val <- tclVar("1")
	rtmz1 <- tkcheckbutton(fr_rtmz, variable = rt.val, command = function(...){
		if(tclvalue(rt.val)=="1"){
			tclvalue(rtmz_choice) = 1
			tkconfigure(rtstart.entry, state = "normal")
			tkconfigure(rtend.entry, state = "normal")
		}else{
			tkconfigure(rtstart.entry, state = "disable")
			tkconfigure(rtend.entry, state = "disable")
		}
	})
	rtstart <- tclVar()
	rtstart.entry <- tkentry(fr_rtmz, width = 5, textvariable = rtstart, bg = "white")
	rtend <- tclVar()
	rtend.entry <- tkentry(fr_rtmz, width = 5, textvariable = rtend, bg = "white")
	
	mz.val <- tclVar("0")
	rtmz2 <- tkcheckbutton(fr_rtmz, variable = mz.val, command = function(...){
		if(tclvalue(mz.val)=="1"){
			tclvalue(rtmz_choice) = 1
			tkconfigure(mzstart.entry, state = "normal")
			tkconfigure(mzend.entry, state = "normal")
		}else{
			tkconfigure(mzstart.entry, state = "disable")
			tkconfigure(mzend.entry, state = "disable")
		}
	})
	
	mzstart <- tclVar()
	mzstart.entry <- tkentry(fr_rtmz, width = 5, textvariable = mzstart, bg = "white")
	mzend <- tclVar()
	mzend.entry <- tkentry(fr_rtmz, width = 5, textvariable = mzend, bg = "white")
	tkconfigure(mzstart.entry, state = "disable")
	tkconfigure(mzend.entry, state = "disable")
	
	tkgrid(tklabel(fr_rtmz, text = "   Retention time (RT) and m/z: "), rtmz_yes, tklabel(fr_rtmz, text = "Yes "), rtmz1, tklabel(fr_rtmz, text = "RT: From"), rtstart.entry, tklabel(fr_rtmz, text = " to"), rtend.entry, tklabel(fr_rtmz, text = "(sec) and "), sticky = "w")
	tkgrid(tklabel(fr_rtmz, text = " "), tklabel(fr_rtmz, text = " "), tklabel(fr_rtmz, text = " "), rtmz2, tklabel(fr_rtmz, text = "m/z: From"), mzstart.entry, tklabel(fr_rtmz, text = " to"), mzend.entry)
	tkgrid(tklabel(fr_rtmz, text = "   "), rtmz_no, tklabel(fr_rtmz, text = "No"))
	tkgrid(fr_rtmz, sticky = "w")
	
	
	fr_filter <- tkframe(dlg)
	peak_filter_choice <- tclVar("2")
	filter_yes <- tkradiobutton(fr_filter, variable = peak_filter_choice, value = "1", command = function(...){
		tkconfigure(peak_filter1, state = "normal")
		if(tclvalue(peak_filter1.val)=="1"){
			tkconfigure(peak_filter1.entry, state = "normal")
		}else{
			tkconfigure(peak_filter1.entry, state = "disable")
		}
		tkconfigure(peak_filter2, state = "normal")
		if(tclvalue(peak_filter2.val)=="1"){
			tkconfigure(peak_filter2.entry, state = "normal")
		}else{
			tkconfigure(peak_filter2.entry, state = "disable")
		}
	})
	filter_no <- tkradiobutton(fr_filter, variable = peak_filter_choice, value = "2", command = function(...){
		tkconfigure(peak_filter1, state = "disable")
		tkconfigure(peak_filter1.entry, state = "disable")
		tkconfigure(peak_filter2, state = "disable")
		tkconfigure(peak_filter2.entry, state = "disable")
	})
	
	peak_filter1.val <- tclVar("0")
	peak_filter1 <- tkcheckbutton(fr_filter, variable = peak_filter1.val, command = function(...){
		if(tclvalue(peak_filter1.val)=="1"){
			tclvalue(peak_filter_choice) = 1
			tkconfigure(peak_filter1.entry, state = "normal")
			fr_temp = tktoplevel(); if(isIcon) tk2ico.set(fr_temp,icon)
			tkwm.title(fr_temp, "Zero")
			tkgrid(tklabel(fr_temp, text = "", height = 0, font = fontIntro_para))
			if(!exists("text_zero")) text_zero <<- tclVar(0)
			zero_lab <- tklabel(fr_temp, text = "   Represent zero: ")
			tk2tip(zero_lab, "To prevent the zero (0) has been transformed to other value\n, like generalized log2 transformation, zero will be 2 and \nfor log transformation, zero will be -Inf.") 
			zero_widget <- tkentry(fr_temp, width = 7, textvariable = text_zero, bg = "white")
			tkgrid(zero_lab, zero_widget, tklabel(fr_temp, text = "    "), sticky = "w")
			tkgrid(tklabel(fr_temp, text = "", height = 0, font = fontIntro_para))
			OK.temp<-tkbutton(fr_temp, text="   OK   ", command=function()tkdestroy(fr_temp))
			tkgrid(OK.temp)
			tkgrid(tklabel(fr_temp, text = "", height = 0, font = fontIntro_para))
		}else{
			tkconfigure(peak_filter1.entry, state = "disable")
		}
	})
	
	peak_filter1.thresh <- tclVar("0.8")
	peak_filter1.entry <- tkentry(fr_filter, width = "4", textvariable = peak_filter1.thresh, bg = "white")
	tkconfigure(peak_filter1.entry, state = "disable")
	peak_filter2.val <- tclVar("0")
	peak_filter2 <- tkcheckbutton(fr_filter, variable = peak_filter2.val, command = function(...){
		if(tclvalue(peak_filter2.val)=="1"){
			tkconfigure(peak_filter2.entry, state = "normal")
			tclvalue(peak_filter_choice) = 1
		}else{
			tkconfigure(peak_filter2.entry, state = "disable")
		}
	})
	
	
	peak_filter2.thresh <- tclVar("0.5")
	peak_filter2.entry <- tkentry(fr_filter, width = "4", textvariable = peak_filter2.thresh, bg = "white")
	tkconfigure(peak_filter2.entry, state = "disable")
	tkgrid(tklabel(fr_filter, text = "   Peak filter: "), filter_yes, tklabel(fr_filter, text = "Yes"), tklabel(fr_filter, text = "("), peak_filter1, tklabel(fr_filter, text = "zero rate >="), peak_filter1.entry, tklabel(fr_filter, text = "(0~1) or "), peak_filter2, tklabel(fr_filter, text = "missing rate >="), peak_filter2.entry, tklabel(fr_filter, text = "(0~1))     "))
	tkgrid(tklabel(fr_filter, text = "   "), filter_no, tklabel(fr_filter, text = "No"))
	tkgrid(fr_filter, sticky = "w")
	tkconfigure(peak_filter1, state = "disable")
	tkconfigure(peak_filter1.entry, state = "disable")
	tkconfigure(peak_filter2, state = "disable")
	tkconfigure(peak_filter2.entry, state = "disable")
	
	tkgrid(tklabel(dlg, text = ""))
	fr_clust <- tkframe(dlg)
	fr_clust.1 <- tkframe(fr_clust)
	fr_clust.2 <- tkframe(fr_clust, relief="ridge", borderwidth=2)
	tkgrid(fr_clust, sticky = "w")
	tkgrid(fr_clust.1, fr_clust.2, sticky = "w")

	clust.label <- tklabel(fr_clust.1, text = "   Clustering:")

	tkgrid(clust.label, sticky = "w")
	tkgrid(tklabel(fr_clust.1, text = ""), sticky = "w")
	tkgrid(tklabel(fr_clust.1, text = ""), sticky = "w")
	tkgrid(tklabel(fr_clust.1, text = ""), sticky = "w")

	fr_metric <- tkframe(fr_clust.2)
	fr_metric.1 <- tkframe(fr_metric)
	fr_metric.2 <- tkframe(fr_metric)
	tkgrid(fr_metric.1, sticky = "w")
	tkgrid(fr_metric.2, sticky = "w")
	metric_p.val <- tclVar("1")
	metric_p <- tkcheckbutton(fr_metric.2, variable = metric_p.val)
	metric_s.val <- tclVar("0")
	metric_s <- tkcheckbutton(fr_metric.2, variable = metric_s.val)
	metric_k.val <- tclVar("0")
	metric_k <- tkcheckbutton(fr_metric.2, variable = metric_k.val)
	tkgrid(tklabel(fr_metric.1, text = "Metric"), sticky = "w")
	tkgrid(tklabel(fr_metric.2, text = "  "), metric_p, tklabel(fr_metric.2, text = "Pearson correlation metricance"), sticky = "w")
	tkgrid(tklabel(fr_metric.2, text = "  "), metric_s, tklabel(fr_metric.2, text = "Spearman correlation metricance"), sticky = "w")
	tkgrid(tklabel(fr_metric.2, text = "  "), metric_k, tklabel(fr_metric.2, text = "Kendall correlation metricance"), sticky = "w")
	
	                                                        
															
	

	fr_method <- tkframe(fr_clust.2)
	fr_method.1 <- tkframe(fr_method)
	fr_method.2 <- tkframe(fr_method)
	tkgrid(fr_method.1, sticky = "w")
	tkgrid(fr_method.2, sticky = "w")
	method_c.val <- tclVar("1")
	method_c <- tkcheckbutton(fr_method.2, variable = method_c.val)
	method_a.val <- tclVar("1")
	method_a <- tkcheckbutton(fr_method.2, variable = method_a.val)
	method_w.val <- tclVar("1")
	method_w <- tkcheckbutton(fr_method.2, variable = method_w.val)
	tkgrid(tklabel(fr_method.1, text = "Linkage: "), sticky = "w")
	tkgrid(tklabel(fr_method.2, text = "  "), method_c, tklabel(fr_method.2, text = "Complete linkage"), sticky = "w")
	tkgrid(tklabel(fr_method.2, text = "  "), method_a, tklabel(fr_method.2, text = "Average linkage"), sticky = "w")
	tkgrid(tklabel(fr_method.2, text = "  "), method_w, tklabel(fr_method.2, text = "Ward's method"), sticky = "w")
	
	                                                            
																
	
	tkgrid(fr_metric, fr_method, sticky = "w")

	tkgrid(tklabel(dlg, text = ""))
	fr_h_method <- tkframe(dlg)
	h_method.val <- tclVar("1")
	h_method_nor <- tkradiobutton(fr_h_method, variable = h_method.val, value = "1")
	h_method_non <- tkradiobutton(fr_h_method, variable = h_method.val, value = "2")
	h_method_boo <- tkradiobutton(fr_h_method, variable = h_method.val, value = "3")
	tkgrid(tklabel(fr_h_method, text = "   s-value thresholding method: "), h_method_nor, tklabel(fr_h_method, text = "Parametric (Normal)"), h_method_non, tklabel(fr_h_method, text = "Nonparametric"), h_method_boo, tklabel(fr_h_method, text = "Bootstrap      "))
	tkgrid(tklabel(fr_h_method, text = "", height = 0, font = fontIntro_para))
	tkgrid(fr_h_method, sticky = "w")
	

	
	
	
	onClust <- function()
	{
		tkdestroy(dlg)
		tkconfigure(tt,cursor="watch")
		
		outpath = paste(tclvalue(textoutput), "/Quality Control", sep = "")
		dir.create(outpath, showWarnings = F)
		
		
		metric = as.numeric(c(tclvalue(metric_p.val), tclvalue(metric_s.val), tclvalue(metric_k.val)))*c(1,2,3); metric <- metric[metric>0]
		method = as.numeric(c(tclvalue(method_c.val), tclvalue(method_a.val), tclvalue(method_w.val)))*c(1,2,3); method <- method[method>0]
		height_method = as.numeric(tclvalue(h_method.val))
		QCtable = read.table(tclvalue(textAbuninput), header = T, quote = "\"", sep = ifelse(grepl(".txt", tclvalue(textAbuninput)), "\t", ","), as.is = T, fill = T, na.string = c("NA", as.character(tclvalue(text_missing))))
		outpath = paste(outpath, "/S", ncol(QCtable)-3, "_P", nrow(QCtable), sep = "")
		dir.create(outpath, showWarnings = F)

		temp = NULL
		if(tclvalue(rtmz_choice) == "1"){
			if(tclvalue(rt.val)=="1"){
				RT_idx = which(QCtable[,3]>=as.numeric(tclvalue(rtstart)) & QCtable[,3]<=as.numeric(tclvalue(rtend)))
				
				temp = c(temp, c(1:nrow(QCtable))[-RT_idx])
			}
			if(tclvalue(mz.val)=="1"){
				MZ_idx = which(QCtable[,2]>=as.numeric(tclvalue(mzstart)) & QCtable[,2]<=as.numeric(tclvalue(mzend)))
				
				temp = c(temp, c(1:nrow(QCtable))[-MZ_idx])
			}
			temp = unique(temp)
		}
		if(length(temp)>0){
			QCrow = nrow(QCtable)
			QCtable = QCtable[-temp,]
			
			
			write.table(QCtable, file = paste(outpath, "/peak_RTMZ_S", ncol(QCtable)-3, "_P", QCrow-length(temp), ".csv", sep = ""), sep = ",", col.names = T, row.names = F, quote = 1)
			tkmessageBox(title = "Quality Control", message = paste("In the region of RT (and MZ), number of peaks is reduced from ", QCrow, " to ", QCrow-length(temp), ".\nThe file is saved as ''peak_RTMZ_S", ncol(QCtable)-3, "_P", QCrow-length(temp), ".csv''.", sep = ""), icon = "info", type = "ok") 
			cat("Quality Control - In the region of RT (and MZ), number of peaks is reduced from ", QCrow, " to ", QCrow-length(temp), ".\nThe file is saved as ''peak_RTMZ_S", ncol(QCtable)-3, "_P", QCrow-length(temp), ".csv''.\n", sep = "") 
		}else{
			
			
			if(tclvalue(rtmz_choice) == "1"){
				tkmessageBox(title = "Quality Control", message = "All peaks are in the region of RT (and MZ).", icon = "info", type = "ok")
				cat("Quality Control - All peaks are in the region of RT (and MZ).\n")
			}
		}

		

			
            ss <- gsub("^\\d{4,8}_(.*?)|(^[A-Z]\\d{4,8})_(.*?)", "\\1", colnames(QCtable)[-(1:3)])
            sample.info = data.frame(Sample = ss, Group = gsub("_[0-9]", "", ss), stringsAsFactors = F)
            ss.rep <- table(sample.info$Group)
            idx <- NULL
            if(min(ss.rep,na.rm=T)==1){
             if(all(ss.rep==1)){
              tkmessageBox(title = "Error", message = "There is only one replicate per sample so that QC procedure CANNOT be executed!", icon = "error", type = "ok")
              cat("Quality Control(Error) - There is only one replicate per sample so that QC procedure CANNOT be executed!\n")
              tkconfigure(tt,cursor="arrow")
              stop()
             } else {
              idx <- which(sample.info$Group%in%names(ss.rep)[ss.rep%in%1])
              QCtable <- QCtable[,-(idx+3)]
             }
            }

			
			samp <- gsub("(.*?)_.*?$","\\1",colnames(QCtable)[-(1:3)])
			rm.peak <- which(apply(QCtable[,-(1:3)],1,function(x) sum(tapply(x,samp,function(y) sum(!is.na(y)))>1))<2) 
			rm.samp <- which(apply(QCtable[,-(1:3)],2,function(x) sum(!is.na(x)))==0) 
			QCtable <- subset(QCtable,subset=!((1:nrow(QCtable))%in%rm.peak),select=setdiff(1:ncol(QCtable),rm.samp+3))

		peak_infor = QCtable[,c(1:3)]

		temp = NULL
		if(tclvalue(peak_filter_choice) == "1"){
			if(tclvalue(peak_filter1.val)=="1"){
				peak_filter1 = apply(QCtable[,-c(1:3)], 1, function(x){
					return(sum(as.numeric(na.omit(x))==as.numeric(tclvalue(text_zero)))/length(as.numeric(na.omit(x))))
				})
				peak_infor = data.frame(peak_infor, zero_rate = peak_filter1)
				temp = c(temp, which(peak_filter1>=as.numeric(tclvalue(peak_filter1.thresh))))
			}
			if(tclvalue(peak_filter2.val)=="1"){
				peak_filter2 = apply(QCtable[,-c(1:3)], 1, function(x){
					return(sum(is.na(x))/length(x))
				})
				peak_infor = data.frame(peak_infor, missing_rate = peak_filter2)
				temp = c(temp, which(peak_filter2>=as.numeric(tclvalue(peak_filter2.thresh))))
			}
			temp = unique(temp)
		}
		if(length(temp)>0){
			QCtable_temp = QCtable[-temp,]
			peak_infor_temp = peak_infor[-temp,]
		}else{
			QCtable_temp = QCtable
			peak_infor_temp = peak_infor
		}
            write.table(QCtable_temp, file = paste(outpath, "/peak_RTMZ_filter_S", ncol(QCtable)-3, "_P", nrow(QCtable)-length(temp), ".csv", sep = ""), sep = ",", col.names = T, row.names = F, quote = 1)
            tkmessageBox(title = "Quality Control", message = sprintf("There are %d samples with two or more replicates and %d peaks in the validated data for subsequent analyses.\n\nThe file is saved as ''peak_RTMZ_filter_S%d_P%d.csv''.",ncol(QCtable)-3,nrow(QCtable)-length(temp),ncol(QCtable)-3,nrow(QCtable)-length(temp)), icon = "info", type = "ok")
            cat("Quality Control - There are ", ncol(QCtable)-3, " samples with two or more replicates and ", nrow(QCtable)-length(temp), " peaks in the validated data for subsequent analyses.\nThe file is saved as ''peak_RTMZ_filter_S", ncol(QCtable)-3, "_P", nrow(QCtable)-length(temp), ".csv''.\n", sep = "")

		rm("QCtable")
		gc(verbose = T)
		
		
		data_QC = QCtable_temp[, -c(1:3)]
		
		ss = gsub("^\\d{4,8}_(.*?)|(^[A-Z]\\d{4,8})_(.*?)", "\\1", colnames(data_QC))
		sample.info = data.frame(Sample = ss, Group = gsub("_[0-9]", "", ss), stringsAsFactors = F)
		num.replicate = max(as.numeric(names(table(table(sample.info$Group)))))
		
		r = r.value(data_QC, num.replicate = num.replicate, index.group = sample.info$Group)
		rm_idx = is.na(r)
		peak_infor_temp = data.frame(peak_infor_temp, r_value = r)
		if(sum(rm_idx)){
			r = na.omit(r)
			data_QC = data_QC[!rm_idx,]
			QCtable_temp = QCtable_temp[!rm_idx,]
			write.table(QCtable_temp, file = paste(outpath, "/peak_RTMZ_filter_rvNA_S", ncol(QCtable)-3, "_P", nrow(QCtable)-length(temp), ".csv", sep = ""), sep = ",", col.names = T, row.names = F, quote = 1)
			
		}
		
		write.table(peak_infor_temp, file = paste(outpath, "/peak_information.csv", sep = ""), sep = ",", col.names = T, row.names = F, quote = 1)
		
		rvalue <- tktoplevel(); if(isIcon) tk2ico.set(rvalue,icon)
		tkwm.title(rvalue, "Please input the cut-off of r-value")
		tkpack(fr_rvalue <- tkframe(rvalue), side = "top")
		replot <- function(...){
			hist(log(r[r>0 & r<10000]), breaks = 50, xaxt = "n", main = "Histogram of r-value (r = SSW / SSB)", xlab = bquote(r-value~(a~scale~of~log[10])), )
			axis(1, at = log(round(exp(seq(log(range(r[r>0 & r<10000])[1]), log(range(r[r>0 & r<10000])[2]), length = 50)), 2)), labels = round(exp(seq(log(range(r[r>0 & r<10000])[1]), log(range(r[r>0 & r<10000])[2]), length = 50)), 2))
		}
		imgrvalue <- tkrplot(fr_rvalue, replot, hscale = 1, vscale = 1)
		tkpack(imgrvalue, side = "top")
		
		tkpack(fr_thresh <- tkframe(rvalue), side = "top", fill = "x")
		tkpack(tklabel(fr_thresh, text = "     r-value cut-off:"), side = "left")
		
		tkpack(fr_thresh_1 <- tkframe(rvalue), side = "top", fill = "x")
		tkpack(fr_thresh_2 <- tkframe(rvalue), side = "top", fill = "x")
		tkpack(fr_thresh_3 <- tkframe(rvalue), side = "top", fill = "x")

		
		temp = unique(c(Inf, 1, round(quantile(r[r<1], probs = seq(0.93, 0.09, by = -0.07), na.rm = T), 2)))
		temp = temp[which(temp!=0)]
		temp = sapply(temp, function(x) ifelse(sum(r<x)<5, "", x))
		for(i in 1:15){
			if(!is.na(temp[i])){
				assign(paste("thresh_", i, sep = ""), tclVar(temp[i]))
			}else{
				assign(paste("thresh_", i, sep = ""), tclVar())
			}
		}
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		label_1 <- tklabel(fr_thresh_1, text = "    1. ", width = 5)
		textthresh_1 <- tkentry(fr_thresh_1, width = "7", textvariable = thresh_1, bg = "white")
		label_2 <- tklabel(fr_thresh_1, text = "2.", width = 3)
		textthresh_2 <- tkentry(fr_thresh_1, width = "7", textvariable = thresh_2, bg = "white")
		label_3 <- tklabel(fr_thresh_1, text = "3.", width = 3)
		textthresh_3 <- tkentry(fr_thresh_1, width = "7", textvariable = thresh_3, bg = "white")
		label_4 <- tklabel(fr_thresh_1, text = "4.", width = 3)
		textthresh_4 <- tkentry(fr_thresh_1, width = "7", textvariable = thresh_4, bg = "white")
		label_5 <- tklabel(fr_thresh_1, text = "5.", width = 3)
		textthresh_5 <- tkentry(fr_thresh_1, width = "7", textvariable = thresh_5, bg = "white")
		tkgrid(label_1, textthresh_1, label_2, textthresh_2, label_3, textthresh_3, label_4, textthresh_4, label_5, textthresh_5, sticky = "w")
		
		label_6 <- tklabel(fr_thresh_2, text = "    6. ", width = 5)
		textthresh_6 <- tkentry(fr_thresh_2, width = "7", textvariable = thresh_6, bg = "white")
		label_7 <- tklabel(fr_thresh_2, text = "7.", width = 3)
		textthresh_7 <- tkentry(fr_thresh_2, width = "7", textvariable = thresh_7, bg = "white")
		label_8 <- tklabel(fr_thresh_2, text = "8.", width = 3)
		textthresh_8 <- tkentry(fr_thresh_2, width = "7", textvariable = thresh_8, bg = "white")
		label_9 <- tklabel(fr_thresh_2, text = "9.", width = 3)
		textthresh_9 <- tkentry(fr_thresh_2, width = "7", textvariable = thresh_9, bg = "white")
		label_10 <- tklabel(fr_thresh_2, text = "10.", width = 3)
		textthresh_10 <- tkentry(fr_thresh_2, width = "7", textvariable = thresh_10, bg = "white")
		tkgrid(label_6, textthresh_6, label_7, textthresh_7, label_8, textthresh_8, label_9, textthresh_9, label_10, textthresh_10, sticky = "w")
		
		label_11 <- tklabel(fr_thresh_3, text = "    11. ", width = 5)
		textthresh_11 <- tkentry(fr_thresh_3, width = "7", textvariable = thresh_11, bg = "white")
		label_12 <- tklabel(fr_thresh_3, text = "12.", width = 3)
		textthresh_12 <- tkentry(fr_thresh_3, width = "7", textvariable = thresh_12, bg = "white")
		label_13 <- tklabel(fr_thresh_3, text = "13.", width = 3)
		textthresh_13 <- tkentry(fr_thresh_3, width = "7", textvariable = thresh_13, bg = "white")
		label_14 <- tklabel(fr_thresh_3, text = "14.", width = 3)
		textthresh_14 <- tkentry(fr_thresh_3, width = "7", textvariable = thresh_14, bg = "white")
		label_15 <- tklabel(fr_thresh_3, text = "15.", width = 3)
		textthresh_15 <- tkentry(fr_thresh_3, width = "7", textvariable = thresh_15, bg = "white")
		tkgrid(label_11, textthresh_11, label_12, textthresh_12, label_13, textthresh_13, label_14, textthresh_14, label_15, textthresh_15, sticky = "w")

		
		tkpack(tklabel(rvalue, text = "", height = 0, font = fontIntro_para), side = "bottom")
		tkpack(fr_conti <- tkframe(rvalue), side = "top", fill = "x")
		tkpack(tkbutton(fr_conti, text = "Continue", width = 20, command = function(...){
			tkdestroy(rvalue)
			tkconfigure(tt,cursor="watch")
			r.range <- sort(unique(as.numeric(na.omit(as.numeric(c(tclvalue(thresh_1), tclvalue(thresh_2), tclvalue(thresh_3), tclvalue(thresh_4), tclvalue(thresh_5), tclvalue(thresh_6), tclvalue(thresh_7), tclvalue(thresh_8), tclvalue(thresh_9), tclvalue(thresh_10), tclvalue(thresh_11), tclvalue(thresh_12), tclvalue(thresh_13), tclvalue(thresh_14), tclvalue(thresh_15)))))), decreasing = T)
			pb <- tkProgressBar("Quality Control - Please wait for QC processing", "0% done", 0, 100, 0, width = 500)
			cat("Quality Control - Please wait for QC processing 0 ")

			outpath.ori <- outpath
			outpath <- sprintf("%s/%s", outpath, c("Parametric", "Nonparametric", "Bootstrap")[height_method])
			for(m in 1:length(metric)){
				acc <- list()
				for(i in 1:length(r.range)){
					acc[[i]] <- r.cluster(data_QC, r = r, cut.r = r.range[i], num.replicate = num.replicate, sample.info.temp = sample.info, c.path = outpath.ori, index.plot = T, index.plot.sh = T, index.correlation = T, index.output = T, height_method = height_method, metric = c("pearson", "spearman", "kendall")[metric[m]], method = c("complete", "average", "ward")[method])

					info <- sprintf("%d%%", round(100*(m/length(metric))*(i/length(r.range))))
					setTkProgressBar(pb, value = round(100*(((m-1)*length(r.range)+i)/(length(metric)*length(r.range)))), sprintf("Quality Control - Please wait for QC processing (%s done)", info), paste(r.range[i], info, sep = " "))
					if(round(100*i/length(r.range))<100){
						cat(round(100*(((m-1)*length(r.range)+i)/(length(metric)*length(r.range)))), " ", sep = "")
					}else{
						cat(round(100*(((m-1)*length(r.range)+i)/(length(metric)*length(r.range)))), " \n", sep = "")
					}
					Sys.sleep(0.1)
				}

				
				if(!all(is.na(sapply(acc, function(x) x$accuracy)))){
					png(sprintf("%s/%s/Accuracy.png", outpath, c("Pearson", "Spearman", "Kendall")[metric[m]]), width = 1440, height = 800)
					
						par(mar = c(4, 5, 4, 2), mgp = c(2.5, 0.6, 0), cex.lab = 2, cex.main = 3, cex.axis = 2, tck = -0.005, font.lab = 2)
						ylim = range(sapply(acc, function(x) x$accuracy),na.rm=T)
						ylim[1] = ylim[1]*0.8
						plot(1:length(acc), type = "n", ylim = c(ylim[1], 1), axes = F, xlab = "r-value", ylab = "Accuracy", main = "Accuracy plot")
						for(i in 1:length(method)){
							acc.temp = sapply(acc, function(x) x$accuracy[[c("complete", "average", "ward")[method[i]]]])
							lines(acc.temp, col = c("red", "blue", "green")[method[i]], pch = c(21, 22, 24)[method[i]], lwd = 2, type = "b", cex = 3)
						}
						abline(h = 1, lty = 2, col = "purple")
						axis(2)
						axis(1, labels = r.range, at = 1:length(acc))
						axis(1, labels = sapply(acc, function(x) x$num.probe), at = 1:length(acc), line = -2, tick = F)
						legend("bottomleft", c("Complete linkage", "Average linkage", "Ward's method")[method], col = c("red", "blue", "green")[method], pch = c(21, 22, 24)[method], lty = 1, horiz = T, lwd = 2, cex = 3, bty = "n", bg = "white") 
					dev.off()
				}
	
				tkconfigure(tt,cursor="arrow")
				outfile <- sprintf("%s/%s/QC_results-%s-%s.RData", outpath, c("Pearson", "Spearman", "Kendall")[metric[m]], Sys.Date(), format(Sys.time(), "%H%M"))
				save(acc, r, r.range, metric, method, m, outpath, height_method, QCtable_temp, file = outfile)
				
				
				

				if((all(is.na(sapply(acc, function(x) x$accuracy))))){
					tkmessageBox(title = "Error", message = "There are no perfectly-clustered samples for the subsequent QC procedure.", icon = "error", type = "ok")
					cat("Quality Control(Error) - There are no perfectly-clustered samples for the subsequent QC procedure.\n")
					stop("Quality Control(Error) - There are no perfectly-clustered samples for the subsequent QC procedure.\n")
				}
				
					
					
					
				
			}
			setTkProgressBar(pb, value = 100, "Quality Control - Please wait for QC processing (100% done)", "Finished 100%")
			Sys.sleep(1)
			close(pb)

			tkmessageBox(title = "Quality Control", message = paste("Output of quality control is finished.\nRData file is saved for Quality Control - Output Viewer", sep = ""), icon = "info", type = "ok")
			
			cat("Quality Control - Output of quality control is finished.\n")

			
				outfiles <- dir(outpath, pattern="\\.RData", full.name=T, recursive=T)
				outfiles.desc <- gsub(".*/(.*?)/(.*?)/(.*?)\\.RData","\\3 (\\1 - \\2CorrDist)",outfiles)

				dlg <- tktoplevel(); if(isIcon) tk2ico.set(dlg,icon)
				tkwm.title(dlg, "QC result")
				listscr <- tkscrollbar(dlg, repeatinterval = 1, command = function(...) tkyview(qc_result, ...))
				qc_result <- tklistbox(dlg, height = max(5,length(outfiles)), selectmode = "single", background = "white", width = 65, yscrollcommand = function(...) tkset(listscr,...))
				sel <<- NULL
				tkbind(qc_result, "<ButtonRelease-1>", function(){
					sel <<- as.integer(tkcurselection(qc_result))+1

					load(outfiles[sel]) 
					QC_interactive(acc = acc, r = r, r.range = r.range, metric = metric[m], method = method, outpath = outpath, height_method = height_method, QCtable_temp = QCtable_temp)
				})
				tkpack(qc_result, side = "left")
				tkpack(listscr, side = "left", fill = "y")
				for(r in 1:length(outfiles)) tkinsert(qc_result, "end", outfiles.desc[r])
			
			
		}), side = "top", anchor = "s")
		if(sum(rm_idx)){
			tkmessageBox(title = "Quality Control", message = paste("The number of peaks with r-value = ''NA'' is ", sum(rm_idx), ".\nThe file is saved as ''peak_RTMZ_filter_rvNA_S", ncol(QCtable)-3, "_P", sum(!rm_idx), ".csv''.", sep = ""), icon = "info", type = "ok")
			cat("Quality Control - The number of peaks with r-value = ''NA'' is ", sum(rm_idx), ".\nThe file is saved as ''peak_RTMZ_filter_rvNA_S", ncol(QCtable)-3, "_P", sum(!rm_idx), ".csv''.\n", sep = "")
		}
		tkmessageBox(title="Quality Control",message="Users can modify the r-value cut-off.", icon = "info", type = "ok")
		cat("Quality Control - Users can modify the r-value cut-off.\n", sep = "")
		tkconfigure(tt,cursor="arrow")

	}
	
	onCancel <- function()
	{
		
		
		tkdestroy(dlg)
		tkfocus(tt)
	}
	tkgrid(tklabel(dlg, text = ""))
	fr <- tkframe(dlg)
	OK.but     <-tkbutton(fr,text="    Run    ",command=onClust)
	Cancel.but <-tkbutton(fr,text="   Cancel   ",command=onCancel)
	tkgrid(tklabel(fr,text="               "), OK.but,tklabel(fr,text="                                      "), Cancel.but, tklabel(fr,text="               "))
	tkgrid(fr)
	tkgrid(tklabel(dlg, text = "", height = 0, font = fontIntro_para))

	tkfocus(dlg)
	tkbind(dlg, "<Destroy>", function() {tkgrab.release(dlg);tkfocus(tt)})
	
	
	tkwait.window(dlg)

	
} 


QualityCtrlViewer <- function()
{
	dlg <- tktoplevel(width = 800); if(isIcon) tk2ico.set(dlg,icon)
	tktitle(dlg) <- "Output Viewer"
	textRDatainput <- tclVar("")
	fr_input <- tkframe(dlg)
	textinputWidget <- tkentry(fr_input,width="50", textvariable = textRDatainput, bg = "white")
	box.input <- tkbutton(fr_input, text = "...",  command = function() tclvalue(textRDatainput) <- tkgetOpenFile(initialfile = tclvalue(textRDatainput), filetypes = "{{RData Files} {.RData}}"))
	
	
	
	
	tkgrid(tklabel(fr_input, text = "", height = 0, font = fontIntro_para))
	tkgrid(tklabel(fr_input, text="   Input .RData file:              "),textinputWidget, box.input, tklabel(fr_input,text="    "), sticky = "w")
	
	tkgrid(tklabel(fr_input, text = "", height = 0, font = fontIntro_para))
	tkgrid(fr_input, sticky = "w")
	
	onOK <- function()
	{
		if(file.exists(tclvalue(textRDatainput))){
			tkdestroy(dlg)
			load(tclvalue(textRDatainput))
			if(all(file.exists(paste(dirname(tclvalue(textRDatainput)), "/", paste(toupper(substring(ls(acc[[1]]$backcol), 1,1)), substring(ls(acc[[1]]$backcol), 2), sep = ""), sep = "")))) outpath = dirname(tclvalue(textRDatainput))

			if(!all(file.exists(paste(outpath, "/", paste(toupper(substring(ls(acc[[1]]$backcol), 1,1)), substring(ls(acc[[1]]$backcol), 2), sep = ""), sep = "")))){
				tkmessageBox(title = "Error", message = paste("The folders of ''", paste(paste(toupper(substring(ls(acc[[1]]$backcol), 1,1)), substring(ls(acc[[1]]$backcol), 2), sep = "")[which(!file.exists(paste(outpath, "/", paste(toupper(substring(ls(acc[[1]]$backcol), 1,1)), substring(ls(acc[[1]]$backcol), 2), sep = ""), sep = "")))], collapse = ", "), "'' not found.\nPlease input the output folder path corresponding to the RData file.", sep = ""), icon = "error", type = "ok")
				
				cat("QC viewer(Redirect Output Folder) - ", paste("The folders of ''", paste(paste(toupper(substring(ls(acc[[1]]$backcol), 1,1)), substring(ls(acc[[1]]$backcol), 2), sep = "")[which(!file.exists(paste(outpath, "/", paste(toupper(substring(ls(acc[[1]]$backcol), 1,1)), substring(ls(acc[[1]]$backcol), 2), sep = ""), sep = "")))], collapse = ", "), " not found. Please input the output folder path corresponding to the RData file.", sep = ""), "\n")
				outpath = tclvalue(tkchooseDirectory())
			}
			mm <- which(sapply(c("Pearson", "Spearman", "Kendall"),regexpr,outpath)>0)
			outpath <- dirname(outpath)
			QC_interactive(acc = acc, r = r, r.range = r.range, metric = metric[mm], method = method, outpath = outpath, height_method = height_method, QCtable_temp = QCtable_temp)
		} else {
			tkmessageBox(title = "Error", message = "RData file is not found.\nPlease input the correct RData file path.", icon = "error", type = "ok")
			cat("QC viewer(Error) - RData file is not found.\nPlease input the correct RData file path.\n")
			tkfocus(dlg)
		}
	}
	
	onCancel <- function()
	{
		
		
		tkdestroy(dlg)
		tkfocus(tt)
	}
	
	fr <- tkframe(dlg)
	OK.but     <-tkbutton(fr,text="     OK     ",command=onOK)
	Cancel.but <-tkbutton(fr,text="   Cancel   ",command=onCancel)
	tkgrid(tklabel(fr,text="    "), OK.but,tklabel(fr,text="           "), Cancel.but, tklabel(fr,text="    "))
	tkgrid(fr)
	tkgrid(tklabel(dlg, text = "", height = 0, font = fontIntro_para))

	tkfocus(dlg)
	tkbind(dlg, "<Destroy>", function() {tkgrab.release(dlg);tkfocus(tt)})
	
	tkwait.window(dlg)

	
}  



DataPrepro <- function()
{
	dlg <- tktoplevel(width = 800); if(isIcon) tk2ico.set(dlg,icon)
	tkwm.title(dlg, "Data Preprocessing")
	
	fr_input <- tkframe(dlg, width = 800)
	
	
	tkgrid(tklabel(fr_input, text = "", height = 0, font = fontIntro_para))
	tkgrid(tklabel(fr_input, text = paste("Input file: ", tclvalue(textAbuninput), "    ", sep = ""), padx = 11, wraplength = 600, justify = "left"), sticky = "w")
	
	
	
	
	tkgrid(tklabel(fr_input, text = paste("Output folder: ", tclvalue(textoutput), "    ", sep = ""), padx = 11, wraplength = 600, justify = "left"), sticky = "w")
	tkgrid(tklabel(fr_input, text = "", height = 0, font = fontIntro_para))
	tkgrid(fr_input, sticky = "w")
	

	ISshift <<- c("None", "Scaling to the mean of IS", "Scaling to the median of IS")
	Trans <<- c("None", "Generalized log2")
	norm.s <<- c("Scale normalization -Mean-", "Scale normalization -Median-", "Quantile normalization", "Standardization", "Pareto scaling")
	norm.p <<- c("Scale normalization -Mean-", "Scale normalization -Median-", "Standardization", "Pareto scaling")
	
	fr_step1 <- tkframe(dlg)
	IS.val = tclVar("None")
	IScomboBox <- ttkcombobox(fr_step1, state = "readonly", values = ISshift, width = 24, textvariable = IS.val)
	tkbind(IScomboBox, "<<ComboboxSelected>>", function(){
		if(tclvalue(IS.val)=="None"){
			tkconfigure(masscombobox, state = "disable")
			tkconfigure(textmass, state = "disable")
		}else{
			tkconfigure(masscombobox, state = "readonly")
			tkconfigure(textmass, state = "normal")
		}
	})
	mass.val = tclVar("m/z")
	masscombobox <- ttkcombobox(fr_step1, state = "readonly", values = c("m/z", "RT"), textvariable = mass.val, width = 6)
	mass <- tclVar("")
	textmass <- tkentry(fr_step1, width = 12, textvariable = mass, bg = "white", validatecommand = "string is double %P", validate = "all")
	tkconfigure(masscombobox, state = "disable")
	tkconfigure(textmass, state = "disable")
	tkgrid(tklabel(fr_step1, text = "   Step 1. Internal standard adjustment:  "), IScomboBox, tklabel(fr_step1, text = "   "), masscombobox, textmass, tklabel(fr_step1, text = "     "))
	tkgrid(fr_step1, sticky = "w")

	Trans.val = tclVar("Generalized log2")
	fr_step2 <- tkframe(dlg)
	TranscomboBox <- ttkcombobox(fr_step2, state = "readonly", values = Trans, width = 24, textvariable = Trans.val)
	tkgrid(tklabel(fr_step2, text = "   Step 2. Transformation:  "), TranscomboBox)
	tkgrid(fr_step2, sticky = "w")
	
	Norm.val <- tclVar("--- select one type of normalizations ---")
	Norm.choose <- tclVar("S") 
	fr_step3 <- tkframe(dlg)
	fr_step3.2 <- tkframe(dlg)
	Norm.rb1.lab <- tklabel(fr_step3,text="None")
	Norm.rb2.lab <- tklabel(fr_step3,text="by Samples")
	Norm.rb3.lab <- tklabel(fr_step3,text="by Peaks")
	Norm.rb1 <- tkradiobutton(fr_step3,command=function() tkconfigure(NormcomboBox,state="disabled"))
	Norm.rb2 <- tkradiobutton(fr_step3,command=function() {tkconfigure(NormcomboBox,value=norm.s,state="normal"); tclvalue(Norm.val) <- "--- select one type of normalizations ---"})
	Norm.rb3 <- tkradiobutton(fr_step3,command=function() {tkconfigure(NormcomboBox,value=norm.p,state="normal"); tclvalue(Norm.val) <- "--- select one type of normalizations ---"})
	NormcomboBox <- ttkcombobox(fr_step3.2, state = "readonly", values = norm.s, width = 32, textvariable = Norm.val)

	tkconfigure(Norm.rb1,variable=Norm.choose,value="N")
	tkconfigure(Norm.rb2,variable=Norm.choose,value="S")
	tkconfigure(Norm.rb3,variable=Norm.choose,value="P")

	tkgrid(tklabel(fr_step3, text = "   Step 3. Normalization:  "),Norm.rb1,Norm.rb1.lab,Norm.rb2,Norm.rb2.lab,Norm.rb3,Norm.rb3.lab)
	tkgrid(tklabel(fr_step3.2, text = "                                      "),NormcomboBox)
	tkgrid.configure(Norm.rb1,Norm.rb1.lab,Norm.rb2,Norm.rb2.lab,Norm.rb3,Norm.rb3.lab,sticky="w")
	tkgrid(fr_step3, sticky = "w")
	tkgrid(fr_step3.2, sticky = "w")

	tkgrid(tklabel(dlg, text = "", height = 1, font = fontIntro_para))
	
	
	onAN <- function()
	{
		temp = TRUE
		if(tclvalue(IS.val) %in% c("1", "2")) temp = ifelse(tclvalue(mass)=="", F, T)
		if(temp){
			ISshiftChoice = which(ISshift == tclvalue(IS.val))-1  
			ISchoice = which(c("m/z", "RT") == tclvalue(mass.val))
			TransChoice = which(Trans == tclvalue(Trans.val))-1  
			
			normalizationDim <- tclvalue(Norm.choose)
			normalizationChoice <- ifelse(normalizationDim=="N",0,which(norm.s == tclvalue(Norm.val)))  


			if(is.na(normalizationChoice)){
				tkmessageBox(title = "Error", message = "Please select one type of normalizations.", icon = "error", type = "ok")
				cat("Data Preprocessing(Error) - Please select one type of normalizations.\n")
				stop("Data Preprocessing(Error) - Please select one type of normalizations.\n")
			}
		
			tkdestroy(dlg)
			tkconfigure(tt,cursor="watch")
			ANtable = read.table(tclvalue(textAbuninput), header = T, sep = ifelse(grepl(".txt", tclvalue(textAbuninput)), "\t", ","), as.is = T, fill = T, na.string = c("NA", as.character(tclvalue(text_missing))), quote = "\"")
			if(tclvalue(mass)!=""){
				ISmass = as.numeric(tclvalue(mass)) 
				ISidx = which.min(abs(ANtable[,ISchoice+1]-ISmass))
				InternalStandard = as.numeric(ANtable[ISidx,-c(1:3)])
				InternalStandard = as.numeric(ANtable[ISidx,-c(1:3)])
				ANtable = ANtable[-ISidx,]
			}
			Peak_Index = ANtable[,1:3]
			ind_name = gsub("^\\d{4,8}_(.*?)|(^[A-Z]\\d{4,8})_(.*?)", "\\1", colnames(ANtable)[-c(1:3)])
			ANdata = ANtable[,-c(1:3)]
			colnames(ANdata) = ind_name

			
			
			if(ISshiftChoice == 1){
				if(any(InternalStandard==0) | any(is.na(InternalStandard))){ 
					tkmessageBox(title = "Error", message = "Zeroes in internal standard.\nAbundance of internal standard should not be zero.", icon = "error", type = "ok")
					cat("Data Preprocessing(Error) - Zeroes in internal standard. Abundance of internal standard should not be zero.\n")
					stop()
				}
				ISadjust = matrix(rep(mean(InternalStandard, na.rm = T)/InternalStandard, nrow(ANdata)), nrow = nrow(ANdata), byrow = T)
				ANdata = ANdata*ISadjust
			}else if(ISshiftChoice == 2){
				if(any(InternalStandard==0) | any(is.na(InternalStandard))){ 
					tkmessageBox(title = "Error", message = "Zeroes in internal standard.\nAbundance of internal standard should not be zero.", icon = "error", type = "ok")
					cat("Data Preprocessing(Error) - Zeroes in internal standard. Abundance of internal standard should not be zero.\n")
					stop()
				}
				ISadjust = matrix(rep(median(InternalStandard, na.rm = T)/InternalStandard, nrow(ANdata)), nrow = nrow(ANdata), byrow = T)
				ANdata = ANdata*ISadjust
			}

			
			
			if(TransChoice == 1){
				ANdata = log2(ANdata+sqrt(ANdata^2+16))
				text_zero <<- tclVar("2")
			}
			
			
			dim <- switch(normalizationDim,N=0,S=2,P=1)
			if(normalizationChoice == 1){
				ANdata = apply(ANdata, dim, function(x) x/mean(x, na.rm = T))
			}else if(normalizationChoice == 2){
				ANdata = apply(ANdata, dim, function(x) x/median(x, na.rm = T))
			}else if(normalizationChoice == 3){
				if(dim==2) ANdata = normalizeQuantiles(as.matrix(ANdata))
				else ANdata = normalizeQuantiles(t(as.matrix(ANdata)))
			}else if(normalizationChoice == 4){
				ANdata = apply(ANdata, dim, function(x) (x-mean(x, na.rm = T))/sd(x, na.rm = T))
			}else if(normalizationChoice == 5){
				ANdata = apply(ANdata, dim, function(x) (x-mean(x, na.rm = T))/sqrt(sd(x, na.rm = T)))
			}
			if(dim==1) ANdata <- t(ANdata)
			
			ANtable = data.frame(Peak_Index, round(ANdata, 5))
			colnames(ANtable)[4:ncol(ANtable)] = ind_name
			outpath = paste(tclvalue(textoutput), "/Data Preprocessing", sep = "")
			dir.create(outpath, showWarnings = F)
			write.table(ANtable, file = paste(outpath, "/", gsub(".*/(.*?)(.txt|.csv)", "\\1", tclvalue(textAbuninput)), "_", c("none", "mean", "median")[ISshiftChoice+1], "_", c("none", "Glog2")[TransChoice+1], "_", c("none", "mean", "median", "QN", "SZ", "PS")[normalizationChoice+1], switch(normalizationDim,N="",S="(S)",P="(P)"), ".csv", sep = ""), col.names = T, row.names = F, quote = 1, sep = ",")
			tkconfigure(tt,cursor="arrow")
			tkmessageBox(title = "Data Preprocessing", message = "Data Preprocessing is done.", icon = "info", type = "ok")
			cat("Data Preprocessing - Data Preprocessing is done.\n", sep = "")
			
			
			
			ReturnVal <- tkmessageBox(title = "Data Preprocessing", message = "Continue to do Quality Control?", icon = "question", type = "yesno", default = "yes")
			if(tclvalue(ReturnVal)=="yes"){
				textAbuninput <<- tclVar(paste(outpath, "/", gsub(".*/(.*?)(.txt|.csv)", "\\1", tclvalue(textAbuninput)), "_", c("none", "mean", "median")[ISshiftChoice+1], "_", c("none", "Glog2")[TransChoice+1], "_", c("none", "mean", "median", "QN", "SZ", "PS")[normalizationChoice+1], switch(normalizationDim,N="",S="(S)",P="(P)"), ".csv", sep = ""))
				QualityControl()
			}
		}else{
			tkmessageBox(title = "Error", message = "M/Z (or RT) should not be empty.\nPlease input the M/Z (or RT) of internal standard.", icon = "error", type = "ok")
			cat("Data Preprocessing(Error) - M/Z (or RT) should not be empty. Please input the M/Z (or RT) of internal standard.\n")
			tkfocus(dlg)
		}
	}
	onCancel <- function()
	{
		
		
		tkdestroy(dlg)
		tkfocus(tt)
	}
	fr <- tkframe(dlg)
	OK.but     <-tkbutton(fr,text="  Run  ",command=onAN)
	Cancel.but <-tkbutton(fr,text="   Cancel   ",command=onCancel)
	tkgrid(tklabel(fr,text="               "), OK.but,tklabel(fr,text="                                            "), Cancel.but, tklabel(fr,text="               "))
	tkgrid(fr)
	tkgrid(tklabel(dlg, text = "", height = 0, font = fontIntro_para))

	tkfocus(dlg)
	tkbind(dlg, "<Destroy>", function() {tkgrab.release(dlg);tkfocus(tt)})
	
	
	tkwait.window(dlg)

	
} 

pcaStat <- function()
{
	dlg <- tktoplevel(width = 800); if(isIcon) tk2ico.set(dlg,icon)
	tkwm.title(dlg, "Batch Effect Detection - PCA")
	
	fr_input <- tkframe(dlg, width = 800)
	
	
		
	
	tkgrid(tklabel(fr_input, text = "", height = 0, font = fontIntro_para), sticky = "we", columnspan = 1)
	tkgrid(tklabel(fr_input, text = paste("Input file: ", tclvalue(textAbuninput), "    ", sep = ""), padx = 11, wraplength = 600, justify = "left"), sticky = "w")
	
	
	
	
	tkgrid(tklabel(fr_input, text = paste("Output folder: ", tclvalue(textoutput), "    ", sep = ""), padx = 11, wraplength = 600, justify = "left"), sticky = "w")
	fr_input.1 <- tkframe(fr_input)
	tkgrid(fr_input.1, sticky = "w")
	if(!exists("textcova")) textcova <<- tclVar("")
	textcovaWidget <- tkentry(fr_input.1,width="50", textvariable = textcova, bg = "white")
	box.cova <- tkbutton(fr_input.1, text = "...",  command = function() tclvalue(textcova) <- tkgetOpenFile(initialfile = as.character(tclvalue(textcova)), filetypes = "{{Text Files} {.txt .csv}}"))
	tkpack(tklabel(fr_input.1,text="   Covariates file:        "), side = "left")
	tkpack(textcovaWidget, side = "left")
	tkpack(box.cova, side = "left")
	tkpack(tklabel(fr_input.1,text="    "), side = "left")
	tkgrid(tklabel(fr_input, text = "", height = 0, font = fontIntro_para), sticky = "w")
	
	tkgrid(fr_input, sticky = "we")
	

	
	onOK <- function()
	{
		if(file.exists(tclvalue(textcova)))
		{
			tkconfigure(dlg,cursor="watch")
			pcatable = read.table(tclvalue(textAbuninput), header = T, fill = T, sep = ifelse(grepl(".txt", tclvalue(textAbuninput)), "\t", ","), na.string = c("", " ", "NA", as.character(tclvalue(text_missing))), quote = "\"")
			Covariate = read.table(tclvalue(textcova), header = T, fill = T, sep = ifelse(grepl(".txt", tclvalue(textcova)), "\t", ","), na.string = c("", " ", "NA", as.character(tclvalue(text_missing))), quote = "\"")
			cova_type = colnames(Covariate)
			cova_type = toupper(substring(cova_type[-1], 1, 1))
			
			Covariate[,1] = gsub("^\\d{4,8}_(.*?)|(^[A-Z]\\d{4,8})_(.*?)", "\\1", Covariate[,1])
			pcadata = pcatable[,-c(1:3)]
			pcanormalize = apply(pcadata, 1, function(x){ x[is.na(x)] = mean(x, na.rm = T);(x-mean(x))/sd(x)})
			colnames(pcadata) = gsub("^\\d{4,8}_(.*?)|(^[A-Z]\\d{4,8})_(.*?)", "\\1", colnames(pcadata))
			tkconfigure(dlg,cursor="arrow")
			if(all(!is.na(match(gsub("_[0-9]", "", colnames(pcadata)), Covariate[,1]))))
			{
				tkdestroy(dlg)
				tkconfigure(tt,cursor="watch")
				mypca = function(pcanormalize){
					pc = list()
					
					
					
					if(nrow(pcanormalize)>=ncol(pcanormalize)){
						pcacor = t(pcanormalize)%*%pcanormalize/nrow(pcanormalize)
						temp1 = eigen(pcacor)
						lambda = temp1$values
						lambda[lambda<0]=0
						
						
						
						pc$scores = pcanormalize%*%temp1$vectors
						pc$scree = lambda
						pc$varprop = pc$scree/sum(pc$scree)
						pc$varpropcum = cumsum(pc$varprop)
					}else{
						pcacor = pcanormalize%*%t(pcanormalize)/nrow(pcanormalize)
						temp1 = eigen(pcacor)
						lambda = temp1$values
						lambda[lambda<0]=0
						
						
						
						pc$scores = temp1$vectors %*% diag(sqrt(lambda*nrow(pcanormalize)))
						pc$scree = lambda
						pc$varprop = pc$scree/sum(pc$scree)
						pc$varpropcum = cumsum(pc$varprop)
					}
					return(pc)
				}
				pcainfo = mypca(pcanormalize)

				Covariate = Covariate[match(gsub("_[0-9]", "", colnames(pcadata)), Covariate[,1]),]
				Covariate = subset(Covariate,select=colnames(Covariate)[-1])
				colnames(Covariate) = substring(colnames(Covariate), 2, nchar(colnames(Covariate)))
				
				PCs <- tktoplevel(); if(isIcon) tk2ico.set(PCs,icon)
				tkwm.title(PCs, "Please input the number of PCs")
				tkpack(fr_pca <- tkframe(PCs), side = "top")
				num_pc <- tclVar("20")
				replot <- function(...){
					par(mfrow = c(1,2))
					plot(pcainfo$scree[1:as.numeric(tclvalue(num_pc))], title = "SCREE plot", type = "l", xlab = "Principal Components", ylab = "Eigenvalue")
					points(pcainfo$scree[1:as.numeric(tclvalue(num_pc))], pch = 21, bg = "lightblue", col = 1)
					abline(h = 1, col = "red")
					par(xpd = NA)
					if(par("usr")[3]<1){
						text(par("usr")[1], 1, label = paste(1, " by KG-rule", sep = ""), pos = 4, col = 3, font = 4, cex = 1)
					}
					text(par("usr")[2], max(pcainfo$scree[1:as.numeric(tclvalue(num_pc))]), label = paste(max(which(pcainfo$scree>1)), " PCs' eigenvalue > 1", sep = ""), pos = 2, col = 3, font = 4, cex = 1)
					
					
					
					plot(pcainfo$varprop[1:as.numeric(tclvalue(num_pc))], title = "Variable Explained", type = "l", xlab = "Principal Components", ylab = "Proportion", ylim = c(0,1))
					points(pcainfo$varprop[1:as.numeric(tclvalue(num_pc))], pch = 21, bg = "red3", col = 1)
					lines(pcainfo$varpropcum[1:as.numeric(tclvalue(num_pc))])
					points(pcainfo$varpropcum[1:as.numeric(tclvalue(num_pc))], pch = 21, bg = "green3", col = 1)
					legend(x = "topright", legend = c("Cumulative", "Proportion"), col = 1, pt.bg = c("green3", "red3"), pch = 21, bty = "o", box.col = "gray", cex = 1)
					text(par("usr")[1], 1, label = paste(sum(pcainfo$varprop>=0.01), " PCs' variance >= 1%", sep = ""), pos = 4, col = 3, font = 4, cex = 1)
					
					
				}
				img <- tkrplot(fr_pca, replot, hscale = 2, vscale = 1.5)
				tkpack(img, side = "top")
				
				tkpack(fr_pc <- tkframe(PCs), side = "top", fill = "x")
				tkpack(tklabel(fr_pc, text = "     The number of the required principal components (PCs):  "), side = "left")
				tkpack(ts_num_pc <- tkscale(fr_pc, width = "14", variable = num_pc, orient = "horizontal", from = 1, to = length(pcainfo$scree), resolution = 1), side = "left")
				tkbind(ts_num_pc, "<ButtonRelease-1>", function(...) tkrreplot(img, hscale = 2, vscale = 1.5))
				tkpack(tklabel(fr_pc, text = " PCs"), side = "left")
				
				tkpack(tklabel(PCs, text = "", height = 0, font = fontIntro_para), side = "bottom")
				tkpack(fr_conti <- tkframe(PCs), side = "top", fill = "x")
				tkpack(tkbutton(fr_conti, text = "  OK  ", width = 10, command = function(...){
					tkdestroy(PCs)
					tkconfigure(tt,cursor="watch")
					outpath = paste(tclvalue(textoutput), "/Batch Effect Detection", sep = "")
					dir.create(outpath, showWarnings = F)
					outpath = paste(outpath, "/Principal Component Analysis (", as.numeric(tclvalue(num_pc)), ")", sep = "")
					dir.create(outpath, showWarnings = F)
					
					png(filename = paste(outpath, "/PCA_screeplot_", tclvalue(num_pc), ".png", sep = ""), width = 1440, height = 800)
						par(mfrow = c(1,2))
						plot(pcainfo$scree[1:as.numeric(tclvalue(num_pc))], main = "Scree Plot", type = "l", xlab = "Principal Components", ylab = "Eigenvalue")
						points(pcainfo$scree[1:as.numeric(tclvalue(num_pc))], pch = 21, bg = "lightblue", col = 1)
						abline(h = 1, col = "red")
						par(xpd = NA)
						if(par("usr")[3]<1){
							text(par("usr")[1], 1, label = paste(1, " by KG-rule", sep = ""), pos = 4, col = 3, font = 4, cex = 1)
						}
						text(par("usr")[2], max(pcainfo$scree[1:as.numeric(tclvalue(num_pc))]), label = paste(max(which(pcainfo$scree>1)), " PCs' eigenvalue > 1", sep = ""), pos = 2, col = 3, font = 4, cex = 1)
						
						
						
						plot(pcainfo$varprop[1:as.numeric(tclvalue(num_pc))], main = "Variation Explained", type = "l", xlab = "Principal Components", ylab = "Proportion", ylim = c(0,1))
						points(pcainfo$varprop[1:as.numeric(tclvalue(num_pc))], pch = 21, bg = "red3", col = 1)
						lines(pcainfo$varpropcum[1:as.numeric(tclvalue(num_pc))])
						points(pcainfo$varpropcum[1:as.numeric(tclvalue(num_pc))], pch = 21, bg = "green3", col = 1)
						legend(x = "topright", legend = c("Cumulative", "Proportion"), col = 1, pt.bg = c("green3", "red3"), pch = 21, bty = "o", box.col = "gray", cex = 1)
						text(par("usr")[1], 1, label = paste(sum(pcainfo$varprop>=0.01), " PCs' variance >= 1%", sep = ""), pos = 4, col = 3, font = 4, cex = 1)
					dev.off()

					pca_scores = as.matrix(pcainfo$scores[,1:as.numeric(tclvalue(num_pc))])
					colnames(pca_scores) = paste("PC_", c(1:as.numeric(tclvalue(num_pc))), sep = "")
					pca_scores = data.frame(SampleID = colnames(pcadata), pca_scores, stringsAsFactors = F)
					pca_var = data.frame(PC_index = 1:as.numeric(tclvalue(num_pc)), PC_eigenvalue = pcainfo$scree[1:as.numeric(tclvalue(num_pc))], PC_varprop = pcainfo$varprop[1:as.numeric(tclvalue(num_pc))], PC_varpropcum = pcainfo$varpropcum[1:as.numeric(tclvalue(num_pc))], stringsAsFactors = F)
					write.table(pca_scores, file = paste(outpath, "/PCA_scores_", as.numeric(tclvalue(num_pc)), ".csv", sep = ""), col.names = T, row.names = F, sep = ",", quote = F)
					write.table(pca_var, file = paste(outpath, "/PCA_variance_", as.numeric(tclvalue(num_pc)), ".csv", sep = ""), col.names = T, row.names = F, sep = ",", quote = F)
					
					
					
					
					for(ind_cov in 1:ncol(Covariate)){
						cova = Covariate[,ind_cov]
						if(cova_type[ind_cov]=="D"){
							color = rainbow(length(unique(as.numeric(factor(na.omit(cova))))))[as.numeric(factor(cova))]
							color[is.na(color)] = "white"
						}else{
							color <- colorRampPalette(c("green3","red"))(length(na.omit(cova)))[rank(cova)]
							color[which(is.na(cova))] = "white"
						}

						png(filename = paste(outpath, "/PCAplot_colourby_", colnames(Covariate)[ind_cov], " (", round(pcainfo$varpropcum[3], 4)*100, "%%).png", sep = ""), width = 1024, height = 768)
							
							if(ind_cov==1){
								
								
								plotpca <<- tktoplevel(); if(isIcon) tk2ico.set(plotpca,icon)
								tkwm.title(plotpca, "PCA plot")
								replot <- function(...){
									if(cova_type[ind_cov]=="D"){
										layout(matrix(1:2,1,),width=c(.9,.1))
										s3d <- scatterplot3d(x = pcainfo$scores[,1], z = pcainfo$scores[,2], y = pcainfo$scores[,3], xlab = paste("Principal Component 1 (", round(pcainfo$varprop[1], 4)*100, "%)", sep = ""), ylab = paste("Principal Component 3 (", round(pcainfo$varprop[3], 4)*100, "%)", sep = ""), zlab = paste("Principal Component 2 (", round(pcainfo$varprop[2], 4)*100, "%)", sep = ""), pch = " ", angle = 50, scale.y = 0.9)

										s3d$points(pcainfo$scores[,1], pcainfo$scores[,3], pcainfo$scores[,2], pch = 21, bg = color)
										par(mar=c(0, 0, 0, 0))
										plot(0,axes=F,type="n")
										legend(x = "right", legend = na.omit(c(levels(factor(na.omit(cova))), ifelse(any(is.na(cova)), "Missing", NA))), pt.bg = na.omit(c(rainbow(length(unique(as.numeric(factor(na.omit(cova)))))), ifelse(any(is.na(cova)), "white", NA))), pch = 21, bty = "o", box.col = "gray", cex = 1.5, col = "black", title=colnames(Covariate)[ind_cov],xpd=NA)
									}else{
										s3d <- scatterplot3d(x = pcainfo$scores[,1], z = pcainfo$scores[,2], y = pcainfo$scores[,3], xlab = paste("Principal Component 1 (", round(pcainfo$varprop[1], 4)*100, "%)", sep = ""), ylab = paste("Principal Component 3 (", round(pcainfo$varprop[3], 4)*100, "%)", sep = ""), zlab = paste("Principal Component 2 (", round(pcainfo$varprop[2], 4)*100, "%)", sep = ""), pch = " ", angle = 50, scale.y = 0.9, mar=c(5, 3, 5, 7)+0.1)

										s3d$points(pcainfo$scores[,1], pcainfo$scores[,3], pcainfo$scores[,2], pch = 21, bg = color)
										par(mar=c(5, 4, 4, 2) + 0.1)
										image.plot(legend.only=TRUE, zlim= c(min(na.omit(cova)), max(na.omit(cova))), nlevel=length(na.omit(cova)), col=colorRampPalette(c("green3","red"))(length(na.omit(cova))),legend.args=list(text = colnames(Covariate)[ind_cov],col="black", cex=2, side=3, line=2))
									}
								}
								imgpca <- tkrplot(plotpca, replot, hscale = 2, vscale = 1.75)
								tkpack(imgpca, side = "top")
								tkpack(tkbutton(plotpca, text = "     OK     ", command = function() tkdestroy(plotpca)), side = "top")
								tkpack(tklabel(plotpca, text = "", height = 0, font = fontIntro_para), side = "bottom")
							}
							if(cova_type[ind_cov]=="D"){
								layout(matrix(1:2,1,),width=c(.9,.1))
								s3d <- scatterplot3d(x = pcainfo$scores[,1], z = pcainfo$scores[,2], y = pcainfo$scores[,3], xlab = paste("Principal Component 1 (", round(pcainfo$varprop[1], 4)*100, "%)", sep = ""), ylab = paste("Principal Component 3 (", round(pcainfo$varprop[3], 4)*100, "%)", sep = ""), zlab = paste("Principal Component 2 (", round(pcainfo$varprop[2], 4)*100, "%)", sep = ""), pch = " ", angle = 50, scale.y = 0.9)

								s3d$points(pcainfo$scores[,1], pcainfo$scores[,3], pcainfo$scores[,2], pch = 21, bg = color)
								par(mar=c(0, 0, 0, .5))
								plot(0,axes=F,type="n")
								legend(x = "right", legend = na.omit(c(levels(factor(na.omit(cova))), ifelse(any(is.na(cova)), "Missing", NA))), col = "black", pt.bg = na.omit(c(rainbow(length(unique(as.numeric(factor(na.omit(cova)))))), ifelse(any(is.na(cova)), "white", NA))), pch = 21, bty = "o", box.col = "gray", cex = switch(ceiling(length(na.omit(c(levels(factor(na.omit(cova))), ifelse(any(is.na(cova)), "Missing", NA))))/25), "1" = 1.5, "2" = 1.2, "3" = 0.9, "4" = 0.6), ncol = ceiling(length(na.omit(c(levels(factor(na.omit(cova))), ifelse(any(is.na(cova)), "Missing", NA))))/25), title=colnames(Covariate)[ind_cov])
							}else{
								s3d <- scatterplot3d(x = pcainfo$scores[,1], z = pcainfo$scores[,2], y = pcainfo$scores[,3], xlab = paste("Principal Component 1 (", round(pcainfo$varprop[1], 4)*100, "%)", sep = ""), ylab = paste("Principal Component 3 (", round(pcainfo$varprop[3], 4)*100, "%)", sep = ""), zlab = paste("Principal Component 2 (", round(pcainfo$varprop[2], 4)*100, "%)", sep = ""), pch = " ", angle = 50, scale.y = 0.9, mar=c(5, 3, 5, 7)+0.1)

								s3d$points(pcainfo$scores[,1], pcainfo$scores[,3], pcainfo$scores[,2], pch = 21, bg = color)
								par(mar=c(5, 4, 4, 2) + 0.1)
								image.plot(legend.only=TRUE, zlim= c(min(na.omit(cova)), max(na.omit(cova))), nlevel=length(na.omit(cova)), col=colorRampPalette(c("green3","red"))(length(na.omit(cova))),legend.args=list(text = colnames(Covariate)[ind_cov],col="black", cex=2, side=3, line=2))
							}
						dev.off()
					}
					tkconfigure(tt,cursor="arrow")
					tkmessageBox(title = "Batch Effect Detection - PCA", message = "Principal Component Analysis is done.", icon = "info", type = "ok")
					cat("Principal Component Analysis-Principal Component Analysis is done.\n")
				}), side = "top", anchor = "s")
				tkfocus(PCs)
				tkconfigure(tt,cursor="arrow")
				tkmessageBox(title="Batch Effect Detection - PCA",message="Please determine the number of principal components (PCs).", icon = "info", type = "ok")
				cat("Principal Components Analysis-Please determine the number of principal components (PCs).\n", sep = "")
			}else{
				tkmessageBox(title = "Error", message = "Sample ID not match.\nPlease input the adequate covariates file.", icon = "error", type = "ok")
				cat("Principal Component Analysis(Error) - Sample ID not match. Please input the adequate covariates file.\n")
				tkfocus(dlg)
			}
		}else{
			tkmessageBox(title = "Error", message = "Covariates file is not found.\nPlease input the correct covariates file path.", icon = "error", type = "ok")
			cat("Principal Component Analysis(Error) - Covariates file is not found. Please input the correct covariates file path.\n")
			tkfocus(dlg)
		}
	}
	onCancel <- function()
	{
		
		
		tkdestroy(dlg)
		tkfocus(tt)
	}
	fr <- tkframe(dlg)
	OK.but     <-tkbutton(fr,text="     OK     ",command=onOK)
	Cancel.but <-tkbutton(fr,text="   Cancel   ",command=onCancel)
	tkgrid(tklabel(fr,text="               "), OK.but,tklabel(fr,text="                                      "), Cancel.but, tklabel(fr,text="               "))
	tkgrid(fr)
	tkgrid(tklabel(dlg, text = "", height = 0, font = fontIntro_para))

	tkfocus(dlg)
	tkbind(dlg, "<Destroy>", function() {tkgrab.release(dlg);tkfocus(tt)})
	
	
	tkwait.window(dlg)

	
}  

LatentGroup <- function()
{
	dlg <- tktoplevel(width = 800); if(isIcon) tk2ico.set(dlg,icon)
	tkwm.title(dlg, "Batch Effect Detection - LG")
	
	fr_input <- tkframe(dlg)
	
	
		
	
	tkgrid(tklabel(fr_input, text = "", height = 0, font = fontIntro_para))
	tkgrid(tklabel(fr_input, text = paste("Input file: ", tclvalue(textAbuninput), "    ", sep = ""), padx = 11, wraplength = 600, justify = "left"), sticky = "w")
	
	
	
	
	tkgrid(tklabel(fr_input, text = paste("Output folder: ", tclvalue(textoutput), "    ", sep = ""), padx = 11, wraplength = 600, justify = "left"), sticky = "w")
		fr_input.1 <- tkframe(fr_input)
	tkgrid(fr_input.1, sticky = "w")
	if(!exists("textcova")) textcova <<- tclVar("")
	textcovaWidget <- tkentry(fr_input.1,width="50", textvariable = textcova, bg = "white")
	box.cova <- tkbutton(fr_input.1, text = "...",  command = function() tclvalue(textcova) <- tkgetOpenFile(initialfile = as.character(tclvalue(textcova)), filetypes = "{{Text Files} {.txt .csv}}"))
	tkpack(tklabel(fr_input.1,text="   Covariates file:        "), side = "left")
	tkpack(textcovaWidget, side = "left")
	tkpack(box.cova, side = "left")
	tkpack(tklabel(fr_input.1,text="    "), side = "left")
	tkgrid(tklabel(fr_input, text = "", height = 0, font = fontIntro_para))
	tkgrid(fr_input, sticky = "w")
	fr_method <- tkframe(dlg)
	method.val <- tclVar("2")
	method_c <- tkradiobutton(fr_method, variable = method.val, value = "1")
	method_a <- tkradiobutton(fr_method, variable = method.val, value = "2")
	method_w <- tkradiobutton(fr_method, variable = method.val, value = "3")
	tkgrid(tklabel(fr_method, text = "   Method:              "), method_c, tklabel(fr_method, text = "Complete linkage"), method_a, tklabel(fr_method, text = "Average linkage"), method_w, tklabel(fr_method, text = "Ward's method"))
	tkgrid(tklabel(fr_method, text = "", height = 0, font = fontIntro_para))
	tkgrid(fr_method, sticky = "w")

	

	
	onOK <- function()
	{
		if(file.exists(tclvalue(textcova))){
			tkconfigure(dlg,cursor="watch")
			BEtable = read.table(tclvalue(textAbuninput), header = T, fill = T, sep = ifelse(grepl(".txt", tclvalue(textAbuninput)), "\t", ","), na.string = c("", " ", "NA", as.character(tclvalue(text_missing))), quote = "\"")
			Covariate = read.table(tclvalue(textcova), header = T, fill = T, sep = ifelse(grepl(".txt", tclvalue(textcova)), "\t", ","), na.string = c("", " ", "NA", as.character(tclvalue(text_missing))), quote = "\"")
			cova_type = colnames(Covariate)
			cova_type = toupper(substring(cova_type[-1], 1, 1))
			
			Covariate[,1] = gsub("^\\d{4,8}_(.*?)|(^[A-Z]\\d{4,8})_(.*?)", "\\1", Covariate[,1])
			outpath = paste(tclvalue(textoutput), "/Batch Effect Detection", sep = "")
			dir.create(outpath, showWarnings = F)
			outpath = paste(outpath, "/Latent Group", sep = "")
			dir.create(outpath, showWarnings = F)
			BEdata = BEtable[,-c(1:3)]
			colnames(BEdata) = gsub("^\\d{4,8}_(.*?)|(^[A-Z]\\d{4,8})_(.*?)", "\\1", colnames(BEdata))
			tkconfigure(dlg,cursor="arrow")
			if(all(!is.na(match(gsub("_[0-9]", "", colnames(BEdata)), Covariate[,1]))))
			{
				tkdestroy(dlg)
				tkconfigure(tt,cursor="watch")
				
				BE_cor = cor(BEdata, method = "pearson", use = "p")
				method = c("complete", "average", "ward")[as.numeric(tclvalue(method.val))]
				hr = hclust(as.dist((1-BE_cor)/2), method)
				method = paste(toupper(substring(method,1,1)), substring(method, 2, nchar(method)), sep = "")
				
				Cutree <- tktoplevel(); if(isIcon) tk2ico.set(Cutree,icon)
				tkwm.title(Cutree, "Batch Effect Detection - LG")
				fr_gap <- tkframe(Cutree)
				png(filename = paste(outpath, "/latent_group_", method, ".png", sep = ""), width = 600, height = 600)
					heatmap(BE_cor, Rowv = as.dendrogram(hr), symm = T, col = colorRampPalette(c("white", "red", "black"), space = "rgb")(20), labRow = NA, labCol = NA, main = "")

					image.plot(legend.only=TRUE, zlim = c(min(na.omit(BE_cor)), max(na.omit(BE_cor))), nlevel=20, col=colorRampPalette(c("white", "red", "black"), space = "rgb")(20), legend.width = 1, legend.mar = 4, axis.args = list(cex.axis = 0.8), legend.args=list(text = "Color\nBar",col="black", cex=0.8, side=3, line=0.5), smallplot=c(.9,.925,.1,.77))
				dev.off()
				
				imgrvalue <<- tclVar()
				tkimage.create("photo", imgrvalue, file = paste(outpath, "/latent_group_", method, ".png", sep = ""))
				tkpack(imgAsLabel <- tklabel(fr_gap, image = imgrvalue))
				tkpack(fr_gap, side = "top")
				tkpack(tklabel(Cutree, text = "", height = 0, font = fontIntro_para), side = "top")
				num_cut = tclVar("2")
				tkpack(fr_cut <- tkframe(Cutree), side = "top", fill = "x")
				tkpack(tklabel(fr_cut, text = "     The number of latent groups:  "), side = "left")
				tkpack(textnumpc <- tkentry(fr_cut, width = "4", textvariable = num_cut, bg = "white", validatecommand = "string is integer %P", xscrollcommand = function(...){
					if(is.numeric(as.numeric(tclvalue(num_cut))) & as.numeric(tclvalue(num_cut))%in%seq.int(2,10)){
						tkconfigure(ok.but, state = "normal")
						k = as.numeric(tclvalue(num_cut))
						cutcl <- cutree(hr, k = k)
						if((k %% 2)==0){
							cutcol_all <- rainbow(2*k-1)
							cutcol = (1+k*(1:k-1)) %% (2*k-1)
							cutcol[cutcol==0] = (2*k-1)
							cutcol = cutcol_all[cutcol]
						}else{
							cutcol_all <- rainbow(k)
							cutcol = (1+(k+1)/2*(1:k-1)) %% k
							cutcol[cutcol==0] = k
							cutcol = cutcol_all[cutcol]
						}
						png(filename = paste(outpath, "/latent_group_", method, "_", k, ".png", sep = ""), width = 600, height = 600)

							heatmap(BE_cor, Rowv = as.dendrogram(hr), symm = T, col = colorRampPalette(c("white", "red", "black"), space = "rgb")(20), labRow = NA, labCol = NA, main = "", ColSideColors = cutcol[cutcl])

							image.plot(legend.only=TRUE, zlim = c(min(na.omit(BE_cor)), max(na.omit(BE_cor))), nlevel=20, col=colorRampPalette(c("white", "red", "black"), space = "rgb")(20), legend.width = 1, legend.mar = 4, axis.args = list(cex.axis = 0.8), legend.args=list(text = "Color\nBar",col="black", cex=0.8, side=3, line=0.5), smallplot=c(.9,.925,.081,.741))
							par(new = T, mar = c(1,1,1,1))
							plot(1:5, type = "n", xlab = "", ylab = "", main = "", axes = F)
							legend(x = "topleft", legend = 1:k, pch = 21, pt.bg = cutcol, border = "white", col = cutcol, ncol = ifelse(k>8, 2, 1), title = "Group index")
						dev.off()
						tkimage.create("photo", imgrvalue, file = paste(outpath, "/latent_group_", method, "_", k, ".png", sep = ""))
						if(exists("label_temp")) tkpack.forget(label_temp)
						
						if(exists("listbox_temp")){tkpack.forget(listbox_temp); tkpack.forget(listscr)}
						if(k>2){
							tkpack(label_temp <<- tklabel(fr_outlier, text = "     Which group indices do you want to remove?  "), side = "left")
							
							listbox_temp <<- tklistbox(fr_outlier, height = 3, selectmode = "multiple", background = "white", width = 5, yscrollcommand = function(...) tkset(listscr,...))
							tkbind(listbox_temp, "<ButtonRelease-1>", function(){
								if(length(as.numeric(tkcurselection(listbox_temp)))>(k-2)){
									tkmessageBox(title = "Error", message = paste("At most remove ", k-2, " groups.\nPlease modify.", sep = ""))
									cat("Latent Group(Error)-", "At most remove ", k-2, " groups.\nPlease modify.\n", sep = "")
									tkconfigure(ok.but, state = "disable")
								}else{
									tkconfigure(ok.but, state = "normal")
								}
							})
							listscr <<- tkscrollbar(fr_outlier, repeatinterval = 1, command = function(...) tkyview(listbox_temp, ...))
							tkpack(listbox_temp, side = "left")
							tkpack(listscr, side = "left", fill = "y")
							for(i in 1:k) tkinsert(listbox_temp, "end", i)
						}
					}else{
						if(as.numeric(tclvalue(num_cut))==""){
							tkmessageBox(title="Error",message="The number of latent groups not found.\nPlease input the number of latent groups.", icon = "error", type = "ok")
							cat("Latent Group(Error)-The number of latent groups not found. Please input the number of sample clusters.\n", sep = "")
						}else{
							tkmessageBox(title="Error",message="The number of latent groups not suitable.\nPlease input the suitable number of latent groups.", icon = "error", type = "ok")
							cat("Latent Group(Error)-The number of latent groups not suitable. Please input the suitable number of latent groups.\n", sep = "")
						}
					}
				}), side = "left")
				tkpack(tklabel(fr_cut, text = " groups (2~10)   "), side = "left")
				
				fr_outlier <<- tkframe(Cutree)
				tkpack(fr_outlier, side = "top", fill = "x")
				tkpack(tklabel(fr_outlier, text = "", height = 0, font = fontIntro_para), side = "top")
				tkpack(tklabel(fr_outlier, text = "", height = 0, font = fontIntro_para), side = "bottom")
				
				tkpack(fr_conti <- tkframe(Cutree), side = "top", fill = "x")
				tkpack(ok.but <<- tkbutton(fr_conti, text = "  OK  ", width = 10, command = function(...){
					tkconfigure(tt,cursor="watch")
					
					k = as.numeric(tclvalue(num_cut))
					if(file.exists(paste(outpath, " (", k, ")", sep = ""))) unlink(paste(outpath, " (", k, ")", sep = ""), recursive = T, force = T)
					file.rename(outpath, paste(outpath, " (", k, ")", sep = ""))
					outpath = paste(outpath, " (", k, ")", sep = "")
					cutcl.ori <- cutcl <- cutree(hr, k = k)
					group.name <- sort(unique(cutcl.ori))
					
					if((k %% 2)==0){
						cutcol_all <- rainbow(2*k-1)
						cutcol = (1+k*(1:k-1)) %% (2*k-1)
						cutcol[cutcol==0] = (2*k-1)
						cutcol = cutcol_all[cutcol]
					}else{
						cutcol_all <- rainbow(k)
						cutcol = (1+(k+1)/2*(1:k-1)) %% k
						cutcol[cutcol==0] = k
						cutcol = cutcol_all[cutcol]
					}
					if(k>2){
						group_out = as.numeric(tkcurselection(listbox_temp))+1
						if(length(group_out)){
							remove_idx = lapply(group_out, function(x) which(cutcl==x))
							cutcl[do.call(c, remove_idx)]=NA
							cutcol[group_out] = "gray60"
							group.name[group_out] <- paste(group.name[group_out],"(NA)")
						}
					}
					tkdestroy(Cutree)
					
					
					cut_output = data.frame(SampleID = names(cutcl), cluster_idx = cutcl, stringsAsFactors = F)
					if(exists("remove_idx")){
						
						write.table(cut_output, file = paste(outpath, "/latent_group", "_", method, "_", k, "_rm", paste(group_out, collapse = "&"), ".csv", sep = ""), col.names = T, row.names = F, sep = ",", quote = F)
					}else{
						write.table(cut_output, file = paste(outpath, "/latent_group", "_", method, "_", k, ".csv", sep = ""), col.names = T, row.names = F, sep = ",", quote = F)
					}

					Covariate = Covariate[match(gsub("_[0-9]", "", cut_output[,1]), Covariate[,1]),]
					Covariate = subset(Covariate,select=colnames(Covariate)[-1])
					colnames(Covariate) = substring(colnames(Covariate), 2, nchar(colnames(Covariate)))
					for(ind_cov in 1:ncol(Covariate)){
						cova = Covariate[,ind_cov]
						if(cova_type[ind_cov]=="D"){
							color <- rainbow(length(unique(as.numeric(factor(na.omit(cova))))))
						}else{
							color <- colorRampPalette(c("green3","red"))(length(na.omit(cova)))
						}
						if(exists("remove_idx")){
							png(filename = paste(outpath, "/latentgroup_colorby_", colnames(Covariate)[ind_cov], "_", method, "_", k, "_rm", paste(group_out, collapse = "&"), ".png", sep = ""), width = 1024, height = 1024)
						}else{
							png(filename = paste(outpath, "/latentgroup_colorby_", colnames(Covariate)[ind_cov], "_", method, "_", k, ".png", sep = ""), width = 1024, height = 1024)
						}
						



							nf <- layout(mat = matrix(c(0,0,2,0,0,3,4,5,1),,3,byrow=TRUE), widths = c(.5,.5,3), height = c(0.25,0.35,3), TRUE)
							
							
							
							
							par(mar = c(3,1,1,3))
							image(1:nrow(BE_cor), 1:nrow(BE_cor), BE_cor[hr$order,hr$order[nrow(BE_cor):1]], col = colorRampPalette(c("white", "red", "black"), space = "rgb")(20),axes = F,xlab="",ylab="",main="",xaxt="n",yaxt="n")
							par(mar = c(0,1,2,3))
							image(1:nrow(BE_cor), c(1,1.3), as.matrix(cutcl.ori[hr$order]),axes = F,xlab="",ylab="",main="",xaxt="n",yaxt="n", col = cutcol)
							par(xpd = NA)
							text(-0.5, 1.15, labels = "latent group", pos = 2, cex = 2)
							par(mar = c(0,1,3.5,3))
							image(1:nrow(BE_cor), c(1,1.3), as.matrix(as.numeric(factor(cova))[hr$order]),axes = F,xlab="",ylab="",xaxt="n",yaxt="n", cex.main = 2, col = color)

							par(xpd = NA)
							text(-0.5, 1.15, labels = colnames(Covariate)[ind_cov], pos = 2, cex = 2)

							par(mar = c(10,5,4.8,6))
							plot(0,type="n",axes=F,main="Group\nindex",cex.main=2,xlab="",ylab="",xlim=c(1,2),ylim=c(40,1),xaxs='i',yaxs='i',xpd=NA)
							sapply(1:length(cutcol),function(x) rect(1,x,1.8,x+1,col=cutcol[x]))
							axis(4,1:length(group.name) + .5,group.name,las=2,tick=F,cex.axis=2)

							if(cova_type[ind_cov]=="D"){
								par(mar = c(10,4,4.8,7))
								plot(0,type="n",axes=F,main=colnames(Covariate)[ind_cov],cex.main=2,xlab="",ylab="",xlim=c(1,2),ylim=c(35,1),xaxs='i',yaxs='i',xpd=NA)
								legend = na.omit(c(levels(factor(na.omit(cova))), ifelse(any(is.na(cova)), "Missing", NA)))
								legend.col = na.omit(c(rainbow(length(unique(as.numeric(factor(na.omit(cova)))))), ifelse(any(is.na(cova)), "white", NA)))
								sapply(1:length(legend),function(x) rect(1,x,1.8,x+1,col=legend.col[x]))
								axis(4,1:length(legend) + .5,legend,las=2,tick=F,cex.axis=2)



							}else{
								par(mar = c(10,4,5,7))
								plot(0,type="n",axes=F,main=colnames(Covariate)[ind_cov],cex.main=2,xlab="",ylab="")
								par(new=T, mar = c(10,0,1,1))
								plot(1:5,type = "n", axes = F,xlab="",ylab="",main="",xaxt="n",yaxt="n")
								image.plot(legend.only=TRUE, zlim = c(min(na.omit(cova)), max(na.omit(cova))), nlevel=length(cova), col=colorRampPalette(c("green3","red"))(length(na.omit(cova))), legend.width = 12, legend.mar = 52, axis.args = list(cex.axis = 2))

							}
						dev.off()
					}
					tkconfigure(tt,cursor="arrow")
					tkmessageBox(title = "Batch Effect Detection - LG", message = "Latent Group analysis is done.", icon = "info", type = "ok")
					cat("Latent Group-Latent Group analysis is done.\n")
				}), side = "top", anchor = "s")
				tkconfigure(ok.but, state = "disable")
				tkpack(tklabel(Cutree, text = "", height = 0, font = fontIntro_para), side = "bottom")
				tkconfigure(tt,cursor="arrow")
				tkmessageBox(title="Batch Effect Detection - LG",message="Please input the number of sample clusters.", icon = "info", type = "ok")
				cat("Latent Group-Please input the number of sample clusters.\n", sep = "")
				tkfocus(Cutree)
			}else{
				tkmessageBox(title = "Error", message = "Sample ID not match.\nPlease input the adequate covariates file.", icon = "error", type = "ok")
				cat("Latent Group(Error) - Sample ID not match. Please input the adequate covariates file.\n")
				tkfocus(dlg)
			}
		}else{
			tkmessageBox(title = "Error", message = "Covariates file is not found.\nPlease input the correct covariates file path.", icon = "error", type = "ok")
			cat("Latent Group(Error) - Covariates file is not found.\nPlease input the correct covariates file path.\n")
			tkfocus(dlg)
		}
	}
	onCancel <- function()
	{
		
		
		tkdestroy(dlg)
		tkfocus(tt)
	}
	fr <- tkframe(dlg)
	OK.but     <-tkbutton(fr,text="     OK     ",command=onOK)
	Cancel.but <-tkbutton(fr,text="   Cancel   ",command=onCancel)
	tkgrid(tklabel(fr,text="               "), OK.but,tklabel(fr,text="                                      "), Cancel.but, tklabel(fr,text="               "))
	tkgrid(fr)
	tkgrid(tklabel(dlg, text = "", height = 0, font = fontIntro_para))

	tkfocus(dlg)
	tkbind(dlg, "<Destroy>", function() {tkgrab.release(dlg);tkfocus(tt)})
	
	
	tkwait.window(dlg)

	
}  


permFact <- function(facts){
	if(NCOL(facts)>1){
		combLevel <- sample(which(!duplicated(facts)))
		idx.combLevel <- tapply(1:nrow(facts),apply(facts,1,paste,collapse="_"),sample)
		
		names(idx.combLevel) <- sample(names(idx.combLevel))
		facts.perm <- facts[rep(combLevel,sapply(idx.combLevel,length))[order(unlist(idx.combLevel),1:nrow(facts))],]
	} else {
		facts.perm <- sample(facts)
	}
	
	return(facts.perm)
}


ANCOVA <- function(){
	dlg <- tktoplevel(); if(isIcon) tk2ico.set(dlg,icon)
	tkwm.title(dlg, "Analysis of Covariance (ANCOVA)")

	fr_input <- tkframe(dlg)
	tkgrid(tklabel(fr_input, text = "", height = 0, font = fontIntro_para))
	tkgrid(tklabel(fr_input, text = paste("Input file: ", tclvalue(textAbuninput), "    ", sep = ""), padx = 11, wraplength = 600, justify = "left"), sticky = "w")
	tkgrid(tklabel(fr_input, text = paste("Output folder: ", tclvalue(textoutput), "    ", sep = ""), padx = 11, wraplength = 600, justify = "left"), sticky = "w")
	
	tkgrid(tklabel(fr_input, text = ""))
	fr_input.1 <- tkframe(fr_input)
	tkgrid(fr_input.1, sticky = "w")
	if(!exists("textCovaInput")) textCovaInput <<- tclVar("")
	textCovaInputWidget <- tkentry(fr_input.1,width="60", textvariable = textCovaInput, bg = "white")
	box.cova <- tkbutton(fr_input.1, text = "...",  command = function() tclvalue(textCovaInput) <- tkgetOpenFile(initialfile = as.character(tclvalue(textCovaInput)), filetypes = "{{Text Files} {.txt .csv}}"))
	cova_label <- tklabel(fr_input.1,text="   Covariates file (optional):       ")
	tk2tip(cova_label, "confounding variable")
	tkgrid(cova_label, textCovaInputWidget, box.cova, tklabel(fr_input.1,text="    "), sticky = "w")
	
	if(!exists("textBatchInput")) textBatchInput <<- tclVar("")
	textbatchWidget <- tkentry(fr_input.1,width="60", textvariable = textBatchInput, bg = "white")
	box.batch <- tkbutton(fr_input.1, text = "...",  command = function() tclvalue(textBatchInput) <- tkgetOpenFile(initialfile = as.character(tclvalue(textBatchInput)), filetypes = "{{Text Files} {.txt .csv}}"))
	batch_label <- tklabel(fr_input.1,text="   Batch effects file (optional):      ")
	tk2tip(batch_label, "from Batch Effect Detection")
	tkgrid(batch_label, textbatchWidget, box.batch, tklabel(fr_input.1,text="    "), sticky = "w")

	if(!exists("textFactorInput")) textFactorInput <<- tclVar("")
	textFactorWidget <- tkentry(fr_input.1,width="60", textvariable = textFactorInput, bg = "white")
	box.factor <- tkbutton(fr_input.1, text = "...",  command = function() tclvalue(textFactorInput) <- tkgetOpenFile(initialfile = as.character(tclvalue(textFactorInput)), filetypes = "{{Text Files} {.txt .csv}}"))
	fac_label <- tklabel(fr_input.1,text="   Factor file (required): ")
	tk2tip(fac_label, "variable of interest")
	tkgrid(fac_label, textFactorWidget, box.factor, tklabel(fr_input.1,text="    "), sticky = "w")
	
	textctrl <- tclVar("")
	tkgrid(tklabel(fr_input.1,text="   Label for the control group:    "), tkentry(fr_input.1,width="9", textvariable = textctrl, bg = "white"), tklabel(fr_input.1,text="    "), sticky = "w")

	tkgrid(tklabel(fr_input, text = ""))
	fr_input.2 <- tkframe(fr_input)
	fr_input.2.1 <- tkframe(fr_input.2)
	fr_input.2.2 <- tkframe(fr_input.2, relief="ridge", borderwidth=2)
	tkgrid(fr_input.2, sticky = "w")
	tkgrid(fr_input.2.1, fr_input.2.2, sticky = "w")

	tkgrid(fr_input.2.1, sticky = "w")
	permute.val <- tclVar(0)
	permute.check <- tkcheckbutton(fr_input.2.1, variable = permute.val, command = function(){
		if(tclvalue(permute.val) == "1"){
			tkconfigure(permutetime.entry, state = "normal")
			tkconfigure(permutecut.entry, state = "normal")
			tkconfigure(permutetime.label, state = "normal")
			tkconfigure(permutecut.label, state = "normal")
			tkconfigure(textcpuWidget, state = "normal")
			tkconfigure(labelcpu, state = "normal")
			tkmessageBox(title = "Warning", message = "Permutation test is time-consuming.", icon = "warning", type = "ok")
			cat("ANCOVA (Warning) - Permutation test is time-consuming.\n")
		}else{
			tkconfigure(permutetime.entry, state = "disable")
			tkconfigure(permutecut.entry, state = "disable")
			tkconfigure(permutetime.label, state = "disable")
			tkconfigure(permutecut.label, state = "disable")
			tkconfigure(textcpuWidget, state = "disabled")
			tkconfigure(labelcpu, state = "disable")
		}
	})

	permute.label <- tklabel(fr_input.2.1, text = "Run permutation test:")

	permutetime.label <- tklabel(fr_input.2.2, text = " Number of permutations:")
	permute.times <- tclVar(10000)
	permutetime.entry <- tkentry(fr_input.2.2, width = 9, textvariable = permute.times, validatecommand = "string is integer %P", validate = "all", bg = "white")
	
	permutecut.label <- tklabel(fr_input.2.2, text = " Cut-off of nominal p-value: ")
	permutecut.val <- tclVar(0.05)
	permutecut.entry <- tkentry(fr_input.2.2, width = 5, textvariable = permutecut.val, validatecommand = "string is double %P", validate = "all", bg = "white")
	
	labelcpu <- tklabel(fr_input.2.2, text = " Number of cores: ")
	tk2tip(labelcpu, "number of cores to be used for permutation test.")
	textcpu <- tclVar(1)
	textcpuWidget <- ttkcombobox(fr_input.2.2, state = "readonly", values = 1:as.numeric(Sys.getenv('NUMBER_OF_PROCESSORS')), width = 4, textvariable = textcpu)

	tkgrid(tklabel(fr_input.2.1,text=" "), permute.check, permute.label, sticky = "w")
	tkgrid(tklabel(fr_input.2.1,text=" "), sticky = "w")
	tkgrid(tklabel(fr_input.2.1,text=" "), sticky = "w")
	tkgrid(tklabel(fr_input.2.2,text=" "), permutetime.label, permutetime.entry, sticky = "w")
	tkgrid(tklabel(fr_input.2.2,text=" "), permutecut.label, permutecut.entry, tklabel(fr_input.2.2,text=" "), sticky = "w")
	tkgrid(tklabel(fr_input.2.2,text=" "), labelcpu, textcpuWidget, sticky  = "w")
	
	tkconfigure(permutetime.entry, state = "disable")
	tkconfigure(permutecut.entry, state = "disable")
	tkconfigure(permutetime.label, state = "disable")
	tkconfigure(permutecut.label, state = "disable")
	tkconfigure(textcpuWidget, state = "disabled")
	tkconfigure(labelcpu, state = "disable")
	
	tkgrid(tklabel(fr_input, text = "", height = 0, font = fontIntro_para))
	tkgrid(fr_input, sticky = "w")

	onOK <- function(){
		check_file_exists <- setNames(c(file.exists(tclvalue(textFactorInput)), tclvalue(textBatchInput)==""|file.exists(tclvalue(textBatchInput)), tclvalue(textCovaInput)==""|file.exists(tclvalue(textCovaInput))), c("Factor", "Batch effect", "Covariates"))
		if(all(check_file_exists)){
			tkconfigure(dlg, cursor = "watch")
			outpath <- paste(tclvalue(textoutput), "/ANCOVA", sep = ""); dir.create(outpath, showWarnings = F)

			ctrl_label <- tclvalue(textctrl)
			ncpu <- as.numeric(tclvalue(textcpu))
			tkconfigure(dlg, cursor = "arrow")

			check_file_exists <- setNames(file.exists(c(tclvalue(textFactorInput), tclvalue(textBatchInput), tclvalue(textCovaInput))), c("factor", "batch effect", "covariates"))
			batch_effect <- Covariate <- NULL

			
			peakabun <- read.table(tclvalue(textAbuninput), header = T, fill = T, sep = ifelse(grepl(".txt", tclvalue(textAbuninput)), "\t", ","), as.is = T, na.string = c("NA", as.character(tclvalue(text_missing))), quote = "\"")  
			peak_index <- peakabun[,c(1:3)]
			peakabun <- peakabun[,-c(1:3)]
			colnames(peakabun) <- gsub("^\\d{4,8}_(.*?)|(^[A-Z]\\d{4,8})_(.*?)", "\\1", colnames(peakabun))

			check_rep <- any(table(gsub("_.*", "", colnames(peakabun)[-(1:3)]))>1) 
			check_id_match <- setNames(rep(TRUE, 3), c("factor", "batch effect", "covariates"))

			
			Factor <- read.table(tclvalue(textFactorInput), header = T, sep = ifelse(grepl(".txt", tclvalue(textFactorInput)), "\t", ","), as.is = T, na.string = c("NA", as.character(tclvalue(text_missing))), quote = "\"", row.names=1)  
			rownames(Factor) <- gsub("^\\d{4,8}_(.*?)|(^[A-Z]\\d{4,8})_(.*?)", "\\1", rownames(Factor))
			
			check_id_match[1] <- !anyNA(match(gsub("_[0-9]", "", colnames(peakabun)), rownames(Factor)))

			
			if(check_file_exists[2]){
				batch_effect <- read.table(tclvalue(textBatchInput), header = T, sep = ifelse(grepl(".txt", tclvalue(textBatchInput)), "\t", ","), as.is = T, na.string = c("NA", as.character(tclvalue(text_missing))), quote = "\"", row.names=1)
				rownames(batch_effect) <- gsub("^\\d{4,8}_(.*?)|(^[A-Z]\\d{4,8})_(.*?)", "\\1", rownames(batch_effect))
				
				check_id_match[2] <- !anyNA(match(colnames(peakabun), rownames(batch_effect)))
			}

			
			if(check_file_exists[3]){
				Covariate <- read.table(tclvalue(textCovaInput), header = T, quote = "\"", fill = T, sep = ifelse(grepl(".txt", tclvalue(textCovaInput)), "\t", ","), as.is = T, na.string = c("NA", as.character(tclvalue(text_missing))),row.names=1)  
				cova_type <- colnames(Covariate)
				cova_type <- toupper(substring(cova_type, 1, 1))
				rownames(Covariate) <- gsub("^\\d{4,8}_(.*?)|(^[A-Z]\\d{4,8})_(.*?)", "\\1", rownames(Covariate))
		
				check_id_match[3] <- !anyNA(match(gsub("_[0-9]", "", colnames(peakabun)), rownames(Covariate)))
			}

			
			if(all(check_id_match)){
				tkdestroy(dlg)
				tkconfigure(tt,cursor="watch")

				peakabun <- t(peakabun)

				DInd <- gsub("_.*", "", rownames(peakabun))
				DRep <- gsub("(.*?)_(.*?)", "\\2", rownames(peakabun))

				Factor <- Factor[match(gsub("_[0-9]", "", rownames(peakabun)), rownames(Factor)), , drop=FALSE]
				if(ctrl_label!=""){ 
					Factor[,1] <- factor(Factor[,1])
					Factor_label <- setdiff(levels(Factor[,1]), ctrl_label)
				} else { 
					Factor_label <- colnames(Factor)
				}

				if(check_file_exists[2]){
					batch_effect <- batch_effect[match(rownames(peakabun), rownames(batch_effect)), , drop=FALSE]
					if(grepl("latent_group_", tclvalue(textBatchInput))) batch_effect[,1] <- factor(batch_effect[,1])
				}

				if(check_file_exists[3]){
					Covariate <- Covariate[match(gsub("_[0-9]", "", rownames(peakabun)), rownames(Covariate)), , drop=FALSE]
					for(x in 1:ncol(Covariate)){
						if(cova_type[x]=="D"){
							Covariate[,x] <- factor(Covariate[,x])
						}
					}
					colnames(Covariate) <- substring(colnames(Covariate), 2, nchar(colnames(Covariate)))
				}

assign("peakabun", peakabun, envir = .GlobalEnv)
assign("Factor", Factor, envir = .GlobalEnv)
assign("batch_effect", batch_effect, envir = .GlobalEnv)
assign("Covariate", Covariate, envir = .GlobalEnv)

				
				
				VIF_table <<- NULL
				vif_next <<- cov_next <<- mod_next <<- FALSE
				if(check_file_exists[2]){
					vif_check <- tktoplevel(); if(isIcon) tk2ico.set(vif_check,icon)
					tkwm.title(vif_check, "Analysis of Covariance - Variance Inflation Factor (VIF)")
					tkgrid(fr_head <- tkframe(vif_check))
					tkgrid(tklabel(fr_head, text = paste("Collinearity Check - VIF", sep = ""), font = fontHeading))
					tkgrid(fr_vif <- tkframe(vif_check))
					
					tclarray <- tclArray()
					VIF_table <- c(ifelse(ctrl_label!="", "Group_comparison", "Quantitative_trait") , "Batch_effect", ifelse(ctrl_label!="", "Group", "Quan"), colnames(Covariate))
					
					for(grp in Factor_label){
						case_idx <- if(ctrl_label!="") which(Factor[,1] %in% c(grp, ctrl_label)) else 1:nrow(Factor)
						
						tmpData <- if(check_file_exists[3]) cbind(Factor[case_idx, 1, drop=FALSE], Covariate[case_idx, , drop=FALSE])
								   else Factor[case_idx, 1, drop=FALSE]
						
						tmp <- t(sapply(1:ncol(batch_effect),function(x){
							sapply(1:ncol(tmpData),function(y){
								tmp <- vif(lm(peakabun[case_idx, 1] ~ tmpData[, y] + batch_effect[case_idx, x]))
								ifelse(class(tmp)=="matrix", round(tmp[2,3], 1), round(tmp[2], 1))
							})
						})) 
						
						VIF_table <- rbind(VIF_table, cbind(ifelse(ctrl_label!="",paste0(grp, "vs.", ctrl_label), grp), colnames(batch_effect), tmp))
					}

					for (i in (1:nrow(VIF_table)))
						for (j in (1:ncol(VIF_table)))
							tclarray[[i-1,j-1]] <- as.tclObj(VIF_table[i,j], drop=T)
					vscr <- tkscrollbar(fr_vif, repeatinterval = 20, orient="vertical", command = function(...) tkyview(table1,...))
					hscr <- tkscrollbar(fr_vif, repeatinterval = 1, orient="horizontal", command = function(...) tkxview(table1,...))
					table1 <- tkwidget(fr_vif,"table",variable=tclarray,rows=nrow(VIF_table),cols=ncol(VIF_table),titlerows=1,titlecols = 2, selectmode="extended",colwidth=20,background="white", yscrollcommand=function(...) tkset(vscr, ...), xscrollcommand = function(...) tkset(hscr, ...))
					tkgrid(table1, vscr); tkgrid.configure(vscr, sticky="nsw")
					tkgrid(hscr, sticky="new")
					tkconfigure(table1,selectmode="extended",rowseparator="\"\n\"",colseparator="\"\t\"")
					tkconfigure(table1,resizeborders="none", state = "disable")

					tkgrid(fr_vif2 <- tkframe(vif_check), sticky = "w")
					tkgrid(tklabel(fr_vif2, text = "", height = 0, font = fontIntro_para))
					vif_threshold <- tclVar("10")
					tkgrid(tklabel(fr_vif2,text="       The VIF upper bound:    "), tkentry(fr_vif2, width="7", textvariable = vif_threshold, bg = "white"), sticky = "w")
				
					tkgrid(fr_next <- tkframe(vif_check))
					tkgrid(tkbutton(fr_next, text = "  Next  ", command = function(...){
						tkdestroy(vif_check)
						tkconfigure(tt,cursor="watch")

						vif_threshold <- as.numeric(tclvalue(vif_threshold))
						write.table(VIF_table, file = paste(outpath, "/VIF_threshold_", vif_threshold, ".csv", sep = ""), sep = ",", col.names = F, row.names = F, quote = F)

						VIF_table <<- matrix(as.numeric(VIF_table[-1,-(1:2)]), nrow = nrow(VIF_table)-1, ncol = ncol(VIF_table)-2, dimnames = list(gsub(paste0("vs.",ctrl_label), "", VIF_table[-1, 1]), VIF_table[1, -(1:2)]))
						VIF_table <<- (VIF_table <= vif_threshold)*1
						
						assign("VIF_table", VIF_table, envir = .GlobalEnv)
						vif_next <<- TRUE
					}))
					tkgrid(tklabel(vif_check, text = ""))
					tkwait.window(vif_check)
				} else {
					vif_next <<- TRUE
				}
				
				
				ancova_cov <<- matrix(0, 1, length(Factor_label), dimnames = list(NULL, Factor_label))
				ancova_cov_idx <<- NULL
				cov_next <<- mod_next <<- FALSE
				if(check_file_exists[3] & vif_next){
					ancova_cov <- outer(colnames(Covariate), Factor_label, function(x, y) sprintf("%s (%s%s)", x, y, ifelse(ctrl_label!="", paste(" vs.",ctrl_label), "")))
					 dimnames(ancova_cov) <- list(colnames(Covariate), Factor_label)
					ancova_var_idx <<- 1:length(ancova_cov)

					cova_select <- tktoplevel(); if(isIcon) tk2ico.set(cova_select,icon)
					tkwm.title(cova_select, "Analysis of Covariance - Covariate Selection")
					choose_frame <- tkframe(cova_select)
					var_frame <- tkframe(choose_frame)
					var_vscr <- tkscrollbar(var_frame, repeatinterval = 5, orient="vertical", command = function(...) tkyview(var_tl,...))
					var_hscr <- tkscrollbar(var_frame, repeatinterval = 1, orient="horizontal", command = function(...) tkxview(var_tl,...))
					var_tl <- tklistbox(var_frame, height = 15, selectmode = "extended", xscrollcommand = function(...) tkset(var_hscr,...), yscrollcommand = function(...) tkset(var_vscr,...), background = "white")
					tkgrid(tklabel(var_frame, text = "Variable"))
					tkgrid(var_tl, var_vscr)
					tkgrid(var_hscr, sticky="new")
					tkgrid.configure(var_vscr, rowspan = 4, sticky = "nsw")
					for(x in ancova_var_idx){
						tkinsert(var_tl, "end", ancova_cov[x])
					}
					cova_frame <- tkframe(choose_frame)
					cova_vscr <- tkscrollbar(cova_frame, repeatinterval = 5, orient="vertical", command = function(...) tkyview(cova_tl,...))
					cova_hscr <- tkscrollbar(cova_frame, repeatinterval = 1, orient="horizontal", command = function(...) tkxview(cova_tl,...))
					cova_tl <- tklistbox(cova_frame, height = 15, selectmode = "extended", xscrollcommand = function(...) tkset(cova_hscr, ...), yscrollcommand = function(...) tkset(cova_vscr, ...), background = "white")
					tkgrid(tklabel(cova_frame, text = "Covariate"))
					tkgrid(cova_tl, cova_vscr)
					tkgrid(cova_hscr, sticky="new")
					tkgrid.configure(cova_vscr, rowspan = 4, sticky = "nsw")
					in_out <- tkframe(choose_frame)
					tkpack(tkbutton(in_out, text = "  >>  ", command = function(...){
						if(length(as.integer(tkcurselection(var_tl)))!=0){
							varIndex <- as.integer(tkcurselection(var_tl))
							for(x in 1:length(varIndex)){
								tkdelete(var_tl, varIndex[x]-x+1)
							}
							ancova_cov_idx <<- sort(c(ancova_cov_idx, ancova_var_idx[varIndex+1]))
							for(x in varIndex){
								tkinsert(cova_tl, which(ancova_cov_idx==ancova_var_idx[x+1])-1, ancova_cov[ancova_var_idx[x+1]])
							}
							ancova_var_idx <<- ancova_var_idx[-(varIndex+1)]
						}
						
					}))
					tkpack(tklabel(in_out, text = "     "))
					tkpack(tkbutton(in_out, text = "  <<  ", command = function(...){
						if(length(as.integer(tkcurselection(cova_tl)))!=0){
							covIndex <- as.integer(tkcurselection(cova_tl))
							for(x in 1:length(covIndex)){
								tkdelete(cova_tl, covIndex[x]-x+1)
							}
							ancova_var_idx <<- sort(c(ancova_var_idx, ancova_cov_idx[covIndex+1]))
							for(x in covIndex){
								tkinsert(var_tl, which(ancova_var_idx==ancova_cov_idx[x+1])-1, ancova_cov[ancova_cov_idx[x+1]])
							}
							ancova_cov_idx <<- ancova_cov_idx[-(covIndex+1)]
						}
						
					}))

					tkgrid(var_frame, in_out, cova_frame, padx = 10)
					tkgrid(tklabel(choose_frame, text = "", height = 0, font = fontIntro_para))
					tkgrid(choose_frame)
					
					tkgrid(fr_next <- tkframe(cova_select))
					tkgrid(but.covSel <- tkbutton(fr_next, text = "  Next  ", state="normal", command = function(...){
					
						tkdestroy(cova_select)
						tkconfigure(tt,cursor="watch")
						if(sum(ancova_cov_idx)){
							ancova_cov[-ancova_cov_idx] <<- 0
							ancova_cov[ancova_cov_idx] <<- 1
						} else {
							ancova_cov[,] <<- 0
						}
						write.table(cbind(c("Cov", rownames(ancova_cov)), rbind(colnames(ancova_cov), ancova_cov)), file = paste(outpath, "/covariate_select.csv", sep = ""), sep = ",", col.names = F, row.names = F, quote = F)

						assign("ancova_cov", ancova_cov, envir = .GlobalEnv)
						cov_next <<- TRUE
					}))
					tkgrid(tklabel(cova_select, text = ""))
					tkwait.window(cova_select)
				} else {
					cov_next <<- TRUE
				}
				
				
				if(vif_next & cov_next){
					mod_next <<- FALSE
					mod.x <<- setNames(vector('list', length(Factor_label)), Factor_label) 
					mod.f <<- setNames(vector('list', length(Factor_label)), Factor_label) 
					for(grp in Factor_label){
						mod.x[[grp]] <<- list(b = NULL, c = NULL)
						if(check_file_exists[3]) mod.x[[grp]][["c"]] <<- rownames(ancova_cov)[ancova_cov[, grp]==1]
						if(check_file_exists[2]) mod.x[[grp]][["b"]] <<- colnames(batch_effect)[apply(VIF_table[rownames(VIF_table)%in%grp, c(ifelse(ctrl_label!="", "Group", "Quan"), mod.x[[grp]][["c"]]), drop = FALSE], 1, function(x) all(x==1))]
					}
					assign("mod.x", mod.x, envir = .GlobalEnv)

					nested_ancova <- tktoplevel(); if(isIcon) tk2ico.set(nested_ancova,icon)
					tkwm.title(nested_ancova, "Analysis of Covariance - Model")
					tkgrid(fr_head <- tkframe(nested_ancova), sticky = "w")
					tkgrid(tklabel(fr_head, text = "Please input an ANCOVA model."))
					tkgrid(fr_input <- tkframe(nested_ancova), sticky = "w")

					tkgrid(tklabel(fr_input, text = "", height = 0, font = fontIntro_para))
					traitname <- ifelse(ctrl_label!="", "Group", "Quan")
					for(grp in Factor_label){
						tkgrid(fr_input.1 <- tkframe(nested_ancova), sticky = "w")
						tkgrid(fr_input.2 <- tkframe(nested_ancova), sticky = "w")
						tkgrid(tklabel(fr_input.1, text = sprintf("  Candidate predictors: %s%s%s%s", traitname, ifelse(check_rep, ", Rep", ""), ifelse(length(mod.x[[grp]][["b"]])!=0, sprintf(", Batch effect (%s)", paste(mod.x[[grp]][["b"]], collapse=", ")), ""), ifelse(length(mod.x[[grp]][["c"]])!=0, sprintf(", Covariates (%s)", paste(mod.x[[grp]][["c"]], collapse=", ")), ""))))
						
						assign(paste0("f_", grp), tkentry(fr_input.2, width = "80", textvariable = tclVar(ifelse(check_rep, paste(traitname, "Rep", sep=":"), "")), bg = "white", state = ifelse(length(mod.x[[grp]][["c"]])!=0 | check_rep, "normal", "disable")))
						tkgrid(tklabel(fr_input.2, text = sprintf("    %s%s%s: %s %s", grp, ifelse(ctrl_label!="", " vs. ", ""), ctrl_label, traitname, ifelse(length(mod.x[[grp]][["b"]])!=0, "+ Batch effect + ", "+ "))), get(paste0("f_", grp)))
					}

					tkgrid(tklabel(nested_ancova, text = "", height = 0, font = fontIntro_para))
					tkgrid(fr_note <- tkframe(nested_ancova), sticky = "w")
					tkgrid(tklabel(fr_note, text = "", height = 0, font = fontIntro_para))
					tkgrid(tklabel(fr_note, text = "     Note: Notation ''A:B'' indicates that variable B nested in variable A.", justify = "left"))
					tkgrid(tklabel(fr_note, text = "     "))
					tkgrid(tkbutton(nested_ancova, text = "  Next  ", command = function(...){
						for(grp in Factor_label){
							mod.f[[grp]] <- gsub("\\+$", "", gsub("\\+{2,}", "\\+", paste(traitname, paste(mod.x[[grp]][["b"]], collapse = "+"), tclvalue(tkget(get(paste0("f_", grp)))), sep = "+")))
						}
						assign("mod.f", mod.f, envir = .GlobalEnv)
						mod_next <<- TRUE

						tkdestroy(nested_ancova)
						tkconfigure(tt,cursor="watch")
					}))
					tkgrid(tklabel(nested_ancova, text = ""))
					tkwait.window(nested_ancova)
				}
				
				if(mod_next){
					prg_all <- 0
					pb <- tkProgressBar("ANCOVA - Please wait for ANCOVA processing", "0% done", 0, 100, 0, width = 500)
					cat("ANCOVA - Please wait for ANCOVA processing 0 ")
					for(grp in Factor_label){
						fa <- which(Factor_label%in%grp)
						case_idx <- if(ctrl_label!="") which(Factor[,1] %in% c(grp, ctrl_label)) else 1:nrow(Factor)
						
						temp_fac <- Factor[case_idx, ifelse(ctrl_label!="", 1, grp)]
						temp_all <- if(check_file_exists[2] & check_file_exists[3]) data.frame(Covariate[case_idx, , drop = FALSE], batch_effect[case_idx, , drop = FALSE], Rep = DRep[case_idx])
									else if (check_file_exists[2] & !check_file_exists[3]) data.frame(batch_effect[case_idx, , drop = FALSE], Rep = DRep[case_idx])
									else if (!check_file_exists[2] & check_file_exists[3]) data.frame(Covariate[case_idx, , drop = FALSE], Rep = DRep[case_idx])
									else data.frame(Rep = DRep[case_idx])

						options(contrasts = c("contr.helmert", "contr.poly"))
						ancova <- apply(peakabun[case_idx,], 2, function(x){
							try(aov(as.formula(sprintf("x ~ %s", gsub("Group|Quan", "temp_fac", mod.f[[grp]]))), data = temp_all), silent=T)
							
						})
						
						ancova.p <- sapply(ancova, function(x){
							tmp <- if(is(x, "aov")) try(Anova(x, type = "III", singular.ok = T)[, "Pr(>F)"], silent = T)
							
							if(!is.numeric(tmp)) tmp <- NA
							
							tmp
						})
						ancova.p <- switch(class(ancova.p), "matrix" = t(ancova.p), "list" = do.call("rbind", ancova.p))
						colnames(ancova.p) <- c("(Intercept)", gsub("temp_fac", sprintf("%s%s%s", grp, ifelse(ctrl_label!="", "_", ""), ctrl_label), colnames(attributes(ancova[[which.min(!is.na(ancova.p[, 1]))]]$terms)$factors)), "")

						options(contrasts = c("contr.SAS", "contr.poly"))
						if(ctrl_label!="") contrasts(temp_fac) = contrasts(temp_fac)[c(which(levels(temp_fac)!=ctrl_label), which(levels(temp_fac)==ctrl_label)),]
						reg <- apply(peakabun[case_idx,], 2, function(x){
							try(lm(as.formula(sprintf("x ~ %s", gsub("Group|Quan", "temp_fac", mod.f[[grp]]))), data = temp_all, na.action = na.exclude), silent=T)
							
							
						})
						prg_grp <- round(100*1/3)
						prg_all <- round(100*(fa-1)/length(Factor_label) + prg_grp/length(Factor_label))
						setTkProgressBar(pb, value = prg_all, sprintf("ANCOVA - Please wait for ANCOVA processing (%d%% done)", prg_all), sprintf("%s%s%s %d%%", grp, ifelse(ctrl_label!="", " - ", ""), ctrl_label, prg_grp))

						
						reg.b <- sapply(reg, function(x){
							if(is(x, "lm")) x$coefficients else NA
						})
						reg.b <- switch(class(reg.b), "matrix" = t(reg.b), "list" = do.call("rbind", reg.b))
						
						reg.r <- sapply(reg, function(x){
							if(is(x, "lm")) resid(x, na.action = na.exclude) else NA
						})
						reg.r <- switch(class(reg.r), "matrix" = t(reg.r), "list" = do.call("rbind", reg.r))
						prg_grp <- round(100*2/3)
						prg_all <- round(100*(fa-1)/length(Factor_label) + prg_grp/length(Factor_label))
						setTkProgressBar(pb, value = prg_all, sprintf("ANCOVA - Please wait for ANCOVA processing (%d%% done)", prg_all), sprintf("%s%s%s %d%%", grp, ifelse(ctrl_label!="", " - ", ""), ctrl_label, prg_grp))

						write.table(cbind(peak_index, ancova.p[, -ncol(ancova.p)]), file = sprintf("%s/pv_%s%s%s_details.csv", outpath, grp, ifelse(ctrl_label!="", " vs. ", ""), ctrl_label), col.names = T, row.names = F, quote = 1, sep = ",")
						write.table(cbind(peak_index, p.adjust(ancova.p[,2], "BH")), file = sprintf("%s/pFDR_%s%s%s.csv", outpath, grp, ifelse(ctrl_label!=""," vs. ",""), ctrl_label), col.names = c(colnames(peak_index), paste0(grp, ifelse(ctrl_label!="","_",""), ctrl_label)), row.names = F, quote = 1, sep = ",")
						write.table(cbind(peak_index, reg.b[,2]), file = sprintf("%s/coef_%s%s%s.csv", outpath, grp, ifelse(ctrl_label!=""," vs. ",""), ctrl_label), col.names = c(colnames(peak_index), paste0(grp, ifelse(ctrl_label!="","_",""), ctrl_label)), row.names = F, quote = 1, sep = ",")
						write.table(cbind(peak_index, reg.r), file = sprintf("%s/resi_%s%s%s_forPLS.csv", outpath, grp, ifelse(ctrl_label!=""," vs. ",""), ctrl_label), col.names = c(colnames(peak_index), rownames(peakabun[case_idx,])), row.names = F, quote = 1, sep = ",")
					
						
						if(tclvalue(permute.val)=="1"){
							peak_sig <- which(ancova.p[, 2]<as.numeric(tclvalue(permutecut.val)))
							if(length(peak_sig)>0){
								pvalue_sig <- ancova.p[peak_sig, 2]
								peakabun_sig <- as.data.frame(peakabun[case_idx, peak_sig])

								cl <- makeCluster(ncpu, type="SOCK")
								clusterExport(cl, c("temp_fac", "peakabun_sig", "pvalue_sig", "temp_all", "mod.f", "grp", "osIsWin"), envir=environment())
								
								registerDoSNOW(cl)
								pvalue_permute <- foreach(i=1:as.numeric(tclvalue(permute.times)), .combine='+', .packages=c("car", "utils", "tcltk")) %dopar% {
								
									pb2 <- tkProgressBar(sprintf("ANCOVA-Parallel task (permutation %d)",i), "0% done", 0, 100, 0, width = 500)
									temp_fac_perm <- sample(temp_fac)
									
										
									
									ancova_perm.p <- apply(peakabun_sig, 2, function(x) {
										tmp <- try(Anova(aov(as.formula(sprintf("x ~ %s", gsub("Group|Quan", "temp_fac_perm", mod.f[[grp]]))), data = temp_all), type = "III", singular.ok = T)[2, "Pr(>F)"], silent=T)

										as.numeric(tmp)
									})
									for(j in 1:100) setTkProgressBar(pb2, value = j, sprintf("ANCOVA-Parallel task (permutation %d)", i), "")
									close(pb2)
									
									ancova_perm.p<pvalue_sig
								}
								stopCluster(cl)

								empvalue <- (pvalue_permute+1)/(as.numeric(tclvalue(permute.times))+1)
								empvalue <- cbind(empvalue, p.adjust(c(empvalue, rep(max(empvalue, as.numeric(tclvalue(permutecut.val))), nrow(peak_index)-length(empvalue))), method = "BH")[1:length(empvalue)])
								write.table(cbind(peak_index[peak_sig,], empvalue), file = sprintf("%s/epv_%s%s%s.csv", outpath, grp, ifelse(ctrl_label!=""," vs. ",""), ctrl_label), col.names = c(colnames(peak_index), paste0(sprintf("%s%s%s",grp,ifelse(ctrl_label!="","_",""),ctrl_label),c("(epv)","(epFDR)"))), row.names = F, quote = 1, sep = ",")
								write.table(cbind(peak_index, reg.r)[peak_sig,][empvalue[, 2]<.05,], file = sprintf("%s/resi_epFDRsig_%s%s%s_forPLS.csv", outpath, grp, ifelse(ctrl_label!=""," vs. ",""), ctrl_label), col.names = c(colnames(peak_index), rownames(peakabun[case_idx,])), row.names = F, quote = 1, sep = ",")
							} else {
								tkmessageBox(title = "Warning", message = "No significant peaks", icon = "warning", type = "ok")
							}
						}

						prg_grp <- round(100*3/3)
						prg_all <- round(100*(fa-1)/length(Factor_label) + prg_grp/length(Factor_label))
						setTkProgressBar(pb, value = prg_all, sprintf("ANCOVA - Please wait for ANCOVA processing (%d%% done)", prg_all), sprintf("%s%s%s %d%%", grp, ifelse(ctrl_label!="", " - ", ""), ctrl_label, prg_grp))
						cat(round(100*fa/length(Factor_label)), ifelse(round(100*fa/length(Factor_label))<100, " ", "\n"), sep = "")
						Sys.sleep(0.1)
					}
					
					setTkProgressBar(pb, value = 100, "ANCOVA - Please wait for ANCOVA processing (100% done)", "Finished 100%")
					Sys.sleep(1)
					close(pb)
					tkconfigure(tt,cursor="arrow")

					toppeak <- tktoplevel(); if(isIcon) tk2ico.set(toppeak,icon)
					tkwm.title(toppeak, "Top peaks")
					fr_top <- tkframe(toppeak)
					tkgrid(fr_top)
					top.lab <- tklabel(fr_top, text = "    How many top peaks do you want to paint red?  ")
					top.num <- tclVar("50")
					top.entry <- tkentry(fr_top, width = 6, textvariable = top.num, bg = "white",validatecommand="string is integer %P", validate = "all")
					
					tkgrid(top.lab, top.entry, tklabel(fr_top, text = paste("(1 ~ ", ncol(peakabun), ")     ", sep = "")))
					tkgrid(tklabel(toppeak, text = "", height = 0, font = fontIntro_para))
					fr_toppeak <- tkframe(toppeak)
					tkgrid(fr_toppeak)
					plot.but <- tkbutton(fr_toppeak, text = "  Plot  ", command = function(...){
						tkconfigure(tt,cursor="watch")
						top.num <- as.numeric(tclvalue(top.num))
						if(!top.num%in%seq.int(1, ncol(peakabun))) top.num = 50
						for(grp in Factor_label){
							pfdr <- read.table(sprintf("%s/pFDR_%s%s%s.csv", outpath, grp, ifelse(ctrl_label!="", " vs. ", ""), ctrl_label), sep = ",", header = T, quote = "\"")
							pvalue <- read.table(sprintf("%s/pv_%s%s%s_details.csv", outpath, grp, ifelse(ctrl_label!="", " vs. ", ""), ctrl_label), sep = ",", header = T, quote = "\"")
							beta_temp <- read.table(sprintf("%s/coef_%s%s%s.csv", outpath, grp, ifelse(ctrl_label!="", " vs. ", ""), ctrl_label), sep = ",", header = T, quote = "\"")

							for(pp in c("pfdr","pvalue")){
								p <- get(pp)
							
								png(filename = sprintf("%s/Volcanoplot_%s%s%s_%s_%d.png", outpath, grp, ifelse(ctrl_label!="", " vs. ", ""), ctrl_label, switch(pp, pfdr="pFDR", pvalue="pv"), top.num), width = 1440, height = 800)
								i <- switch(pp, pfdr=4, pvalue=5)
								par(mar = c(5, 4, 1, 3))
								pch <- rep(21, nrow(p))
								pch[which(-log10(p[,i])>=-log10(0.05))] <- 24
								bg <- rep("white", nrow(p))
								bg[which(-log10(p[,i])>=-log10(0.05))] <- "green3"
								bg[order(-log10(p[,i]), decreasing = T)[1:top.num]] <- "red"
								plot(beta_temp[, 4], -log10(p[,i]), xlab = "", ylab = "", cex = 1, pch = pch, bg = bg)
								text(beta_temp[order(-log10(p[,i]), decreasing = T)[1:top.num], 4], -log10(p[order(-log10(p[,i]), decreasing = T)[1:top.num],i]), label = peak_index[order(-log10(p[,i]), decreasing = T)[1:top.num],1], cex = 1, pos = 1, adj = 0.2)
								abline(h = -log10(0.05), col = "green3")
								par(xpd = NA)
								text(par("usr")[1], -log10(0.05), col = "green3", label = substitute(paste(plain(-log[10]), "(0.05) (", test, " significant peaks)"), list(test = sum(p[,i]<0.05,na.rm=T))), font = 4, cex = 1.5, adj = c(0,-0.1))
								par(xpd = F)
								mtext("Fold Change (coefficient)", line = 2, side = 1, col="black", padj = 1, cex = 1.5)
								mtext(bquote(-log[10](.(switch(pp, pfdr="pFDR", pvalue="pv")))), line = 3.5, side = 2, col="black", padj = 1, cex = 1.5)
								legend(x = "topleft", legend = c(sprintf("%s<0.05", switch(pp, pfdr="pFDR", pvalue="pv")), paste0("top", top.num)), pch = 24, bty = "o", box.col = "black", cex = 1, pt.bg = c("green3", "red"), bg = rgb(1, 1, 1, alpha = 0.7))
								dev.off()
							}
							
							if(tclvalue(permute.val)=="1" & file.exists(sprintf("%s/epv_%s%s%s.csv", outpath, grp, ifelse(ctrl_label!="", " vs. ", ""), ctrl_label))){
								beta_temp <- beta_temp[which(pvalue[, 5]<as.numeric(tclvalue(permutecut.val))),]
								peak_index_temp <- peak_index[which(pvalue[, 5]<as.numeric(tclvalue(permutecut.val))),]
								p <- read.table(sprintf("%s/epv_%s%s%s.csv", outpath, grp, ifelse(ctrl_label!="", " vs. ", ""), ctrl_label), sep = ",", header = T, quote = "\"")

								for(pp in c("pfdr","pvalue")){
									png(filename = sprintf("%s/Volcanoplot_%s%s%s_%s_%d.png", outpath, grp, ifelse(ctrl_label!="", " vs. ", ""), ctrl_label, switch(pp, pfdr="epFDR", pvalue="epv"), top.num), width = 1440, height = 800)
									i <- switch(pp, pfdr=5, pvalue=4)
									par(mar = c(5, 4, 1, 3))
									pch <- rep(21, nrow(p))
									pch[which(-log10(p[,i])>=-log10(0.05))] <- 24
									bg <- rep("white", nrow(p))
									bg[which(-log10(p[,i])>=-log10(0.05))] <- "green3"
									bg[order(-log10(p[,i]), decreasing = T)[1:min(top.num, nrow(p))]] <- "red"
									plot(beta_temp[, 4], -log10(p[,i]), xlab = "", ylab = "", cex = 1, pch = pch, bg = bg)
									text(beta_temp[order(-log10(p[,i]), decreasing = T)[1:min(top.num, nrow(p))], 4], -log10(p[order(-log10(p[,i]), decreasing = T)[1:min(top.num, nrow(p))], i]), label = peak_index_temp[order(-log10(p[,i]), decreasing = T)[1:min(top.num, nrow(p))], 1], cex = 1, pos = 1, adj = 0.2)
									abline(h = -log10(0.05), col = "green3")
									par(xpd = NA)
									text(par("usr")[1], -log10(0.05), col = "green3", label = substitute(paste(plain(-log[10]), "(0.05) (", test, " significant peaks)"), list(test = sum(p[,i]<0.05,na.rm=T))), font = 4, cex = 1.5, adj = c(0,-0.1))
									par(xpd = F)
									mtext("Fold Change (coefficient)", line = 2, side = 1, col="black", padj = 1, cex = 1.5)
									mtext(bquote(-log[10](.(switch(pp, pfdr="epFDR", pvalue="epv")))), line = 3.5, side = 2, col="black", padj = 1, cex = 1.5)
									legend(x = "topleft", legend = c(sprintf("%s<0.05", switch(pp, pfdr="pFDR", pvalue="pv")), paste0("top", min(top.num, nrow(p)))), pch = 24, bty = "o", box.col = "black", cex = 1, pt.bg = c("green3", "red"), bg = rgb(1, 1, 1, alpha = 0.7))
									dev.off()
								}
							}
						}
						tkmessageBox(title = "Analysis of Covariance", message = "Graphs are saved.", icon = "info", type = "ok")
						cat("ANCOVA - Graphs are saved.\n", sep = "")
						tkconfigure(tt,cursor="arrow")
					})
					close.but <- tkbutton(fr_toppeak, text = "  Close  ", command = function(...){
						tkdestroy(toppeak)
						tkmessageBox(title="Analysis of Covariance", message = "ANCOVA is done.", icon = "info", type = "ok")
						cat("ANCOVA - ANCOVA is done.\n", sep = "")
					})
					tkgrid(plot.but, tklabel(fr_toppeak, text = "                              "), close.but)
					tkgrid(tklabel(toppeak, text = "", height = 0, font = fontIntro_para))
				}
			} else {
				tkmessageBox(title = "Error", message = sprintf("Sample IDs are unmatched.\nPlease input the adequate file%s (%s).",ifelse(sum(!check_id_match)>1,"s",""),paste(names(check_id_match)[!check_id_match], collapse = ", ")), icon = "error", type = "ok")
				cat("ANCOVA(Error) - ", sprintf("Sample IDs are unmatched. Please input the adequate file%s (%s).",ifelse(sum(!check_id_match)>1,"s",""),paste(names(check_id_match)[!check_id_match], collapse = ", ")), "\n")
				tkfocus(dlg)
			}
		} else {
			tkmessageBox(title = "Error", message = sprintf("Input file%s (%s) %s not found.\nPlease input the correct file path%s.",ifelse(sum(!check_file_exists)>1,"s",""),paste(names(check_file_exists)[!check_file_exists],collapse=", "),ifelse(sum(!check_file_exists)>1,"are","is"),ifelse(sum(!check_file_exists)>1,"s","")), icon = "error", type = "ok")
			cat("ANCOVA(Error) - ", sprintf("Input file%s (%s) %s not found. Please input the correct file path%s.",ifelse(sum(!check_file_exists)>1,"s",""),paste(names(check_file_exists)[!check_file_exists],collapse=", "),ifelse(sum(!check_file_exists)>1,"are","is"),ifelse(sum(!check_file_exists)>1,"s","")), "\n")
			tkfocus(dlg)
		}
	}

	onCancel <- function(){
		tkdestroy(dlg)
		tkfocus(tt)
	}
	
	tkgrid(tklabel(dlg, text = ""))
	fr <- tkframe(dlg)
	OK.but     <- tkbutton(fr,text="     OK     ",command=onOK,state="normal")
	Cancel.but <- tkbutton(fr,text="   Cancel   ",command=onCancel)
	tkgrid(tklabel(fr,text="               "), OK.but,tklabel(fr,text="                                      "), Cancel.but, tklabel(fr,text="               "))
	tkgrid(fr)
	tkgrid(tklabel(dlg, text = "", height = 0, font = fontIntro_para))

	tkfocus(dlg)
	tkbind(dlg, "<Destroy>", function() {tkgrab.release(dlg);tkfocus(tt)})
	tkwait.window(dlg)
}

ANCOVA_old <- function()
{
	dlg <- tktoplevel(); if(isIcon) tk2ico.set(dlg,icon)
	tkwm.title(dlg, "Analysis of Covariance (ANCOVA)")
	
	fr_input <- tkframe(dlg)
	
	
		
		
	
	tkgrid(tklabel(fr_input, text = "", height = 0, font = fontIntro_para))
	tkgrid(tklabel(fr_input, text = paste("Input file: ", tclvalue(textAbuninput), "    ", sep = ""), padx = 11, wraplength = 600, justify = "left"), sticky = "w")
	
	
	
	
	tkgrid(tklabel(fr_input, text = paste("Output folder: ", tclvalue(textoutput), "    ", sep = ""), padx = 11, wraplength = 600, justify = "left"), sticky = "w")
	
	fr_input.1 <- tkframe(fr_input)
	tkgrid(fr_input.1, sticky = "w")
	if(!exists("textcova")) textcova <<- tclVar("")
	textcovaWidget <- tkentry(fr_input.1,width="60", textvariable = textcova, bg = "white")
	box.cova <- tkbutton(fr_input.1, text = "...",  command = function() tclvalue(textcova) <- tkgetOpenFile(initialfile = as.character(tclvalue(textcova)), filetypes = "{{Text Files} {.txt .csv}}"))
	cov_label <- tklabel(fr_input.1,text="   Covariates file (optional):       ")
	tk2tip(cov_label, "confounding variable")
	tkgrid(cov_label, textcovaWidget, box.cova, tklabel(fr_input.1,text="    "), sticky = "w")
	
	if(!exists("textbatchinput")) textbatchinput <<- tclVar("")
	textbatchWidget <- tkentry(fr_input.1,width="60", textvariable = textbatchinput, bg = "white")
	box.batch <- tkbutton(fr_input.1, text = "...",  command = function() tclvalue(textbatchinput) <- tkgetOpenFile(initialfile = as.character(tclvalue(textbatchinput)), filetypes = "{{Text Files} {.txt .csv}}"))
	batch_label <- tklabel(fr_input.1,text="   Batch effects file (optional):      ")
	tk2tip(batch_label, "from Batch Effect Detection")
	tkgrid(batch_label, textbatchWidget, box.batch, tklabel(fr_input.1,text="    "), sticky = "w")

	if(!exists("textFactorinput")) textFactorinput <<- tclVar("")
	textFactorWidget <- tkentry(fr_input.1,width="60", textvariable = textFactorinput, bg = "white")
	box.factor <- tkbutton(fr_input.1, text = "...",  command = function() tclvalue(textFactorinput) <- tkgetOpenFile(initialfile = as.character(tclvalue(textFactorinput)), filetypes = "{{Text Files} {.txt .csv}}"))
	fac_label <- tklabel(fr_input.1,text="   Factor file (required): ")
	tk2tip(fac_label, "variable of interest")
	tkgrid(fac_label, textFactorWidget, box.factor, tklabel(fr_input.1,text="    "), sticky = "w")
	
	fr_input.2 <- tkframe(fr_input)
	tkgrid(fr_input.2, sticky = "w")
	textctrl <- tclVar("")
	tkgrid(tklabel(fr_input.2,text="   Label for the control group:    "), tkentry(fr_input.2,width="9", textvariable = textctrl, bg = "white"), tklabel(fr_input.2,text="    "), sticky = "w")

	
	fr_input.3 <- tkframe(fr_input)
	tkgrid(fr_input.3, sticky = "w")
	permute.val <- tclVar(0)
	permute.check <- tkcheckbutton(fr_input.3, variable = permute.val, command = function(){
		if(tclvalue(permute.val) == "1"){
			tkconfigure(permute.entry, state = "normal")
			tkconfigure(permutecut.entry, state = "normal")
			tkconfigure(permute.label, state = "normal")
			tkconfigure(permutecut.label, state = "normal")
			tkconfigure(textcpuWidget, state = "normal")
			tkconfigure(labelcpu, state = "normal")
			tkmessageBox(title = "Warning", message = "Permutation test is time-consuming.", icon = "warning", type = "ok")
			cat("ANCOVA (Warning) - Permutation test is time-consuming.\n")
		}else{
			tkconfigure(permute.entry, state = "disable")
			tkconfigure(permutecut.entry, state = "disable")
			tkconfigure(permute.label, state = "disable")
			tkconfigure(permutecut.label, state = "disable")
			tkconfigure(textcpuWidget, state = "disabled")
			tkconfigure(labelcpu, state = "disable")
		}
	})
	permute.label <- tklabel(fr_input.3, text = "Run permutation test      num. of permutations: ")
	permute.times <- tclVar(10000)
	permute.entry <- tkentry(fr_input.3, width = 9, textvariable = permute.times, validatecommand = "string is integer %P", validate = "all", bg = "white")
	
	permutecut.label <- tklabel(fr_input.3, text = "Cut-off of nominal p-value: ")
	permutecut.val <- tclVar(0.05)
	permutecut.entry <- tkentry(fr_input.3, width = 5, textvariable = permutecut.val, validatecommand = "string is double %P", validate = "all", bg = "white")
	tkgrid(tklabel(fr_input.3,text=" "), permute.check, permute.label, permute.entry, permutecut.label, permutecut.entry, tklabel(fr_input.3, text = "    "), sticky = "w")
	
	fr_input.4 <- tkframe(fr_input)
	tkgrid(fr_input.4, sticky = "w")
	labelcpu <- tklabel(fr_input.4, text = "Number of cores: ")
	tk2tip(labelcpu, "number of cores to be used for permutation test.")
	textcpu <- tclVar(1)
	textcpuWidget <- ttkcombobox(fr_input.4, state = "readonly", values = 1:as.numeric(Sys.getenv('NUMBER_OF_PROCESSORS')), width = 4, textvariable = textcpu)
	tkgrid(tklabel(fr_input.4,text=" "), labelcpu, textcpuWidget, sticky  = "w")
	
	tkconfigure(permute.entry, state = "disable")
	tkconfigure(permutecut.entry, state = "disable")
	tkconfigure(permute.label, state = "disable")
	tkconfigure(permutecut.label, state = "disable")
	tkconfigure(textcpuWidget, state = "disabled")
	tkconfigure(labelcpu, state = "disable")
	
	tkgrid(tklabel(fr_input, text = "", height = 0, font = fontIntro_para))
	tkgrid(fr_input, sticky = "w")
	
	
	
	onOK <- function()
	{
		if(any(tclvalue(textcova)=="", all(!(tclvalue(textcova)==""), file.exists(tclvalue(textcova)))) & file.exists(tclvalue(textFactorinput)) & any(tclvalue(textbatchinput)=="", all(!(tclvalue(textbatchinput)==""), file.exists(tclvalue(textbatchinput))))){
			tkconfigure(dlg, cursor = "watch")
			outpath = paste(tclvalue(textoutput), "/ANCOVA", sep = "")
			dir.create(outpath, showWarnings = F)
			peakabun = read.table(tclvalue(textAbuninput), header = T, fill = T, sep = ifelse(grepl(".txt", tclvalue(textAbuninput)), "\t", ","), as.is = T, na.string = c("NA", as.character(tclvalue(text_missing))), quote = "\"")  
			Factor = read.table(tclvalue(textFactorinput), header = T, sep = ifelse(grepl(".txt", tclvalue(textFactorinput)), "\t", ","), as.is = T, na.string = c("NA", as.character(tclvalue(text_missing))), quote = "\"")  
			ctrl_label = tclvalue(textctrl)
			ncpu = as.numeric(tclvalue(textcpu))
			tkconfigure(dlg, cursor = "arrow")

			check <- all(table(gsub("_.*", "", colnames(peakabun)[-(1:3)]))==1)
			
			if((tclvalue(textbatchinput)=="") & (tclvalue(textcova)=="")){  
				peak_index = peakabun[,c(1:3)]
				peakabun = peakabun[,-c(1:3)]
				colnames(peakabun) = gsub("^\\d{4,8}_(.*?)|(^[A-Z]\\d{4,8})_(.*?)", "\\1", colnames(peakabun))
				Factor[,1] = gsub("^\\d{4,8}_(.*?)|(^[A-Z]\\d{4,8})_(.*?)", "\\1", Factor[,1])
				
				if(all(!is.na(match(gsub("_[0-9]", "", colnames(peakabun)), Factor[,1]))))
				{
					tkdestroy(dlg)
				
					tkconfigure(tt,cursor="watch")
					
					Rep = gsub("(.*?)_(.*?)", "\\2", colnames(peakabun))
					Ind <- gsub("_.*", "", colnames(peakabun))
					if(ctrl_label!=""){
						Factor = Factor[match(gsub("_[0-9]", "", colnames(peakabun)), Factor[,1]),2]
						Factor_label = levels(factor(Factor))
						Factor_label = Factor_label[-which(Factor_label==ctrl_label)]
					} else {
						Factor <- Factor[match(gsub("_[0-9]", "", colnames(peakabun)), Factor[,1]),]
						Factor <- subset(Factor,select=-1) 
						Factor_label <- colnames(Factor)
					}

					peakabun = t(peakabun)

					
					op <- options(contrasts = c("contr.helmert", "contr.poly"))
					pb <- tkProgressBar("ANCOVA-Please wait for ANCOVA processing", "0% done", 0, 100, 0, width = 500)
					cat("ANCOVA-Please wait for ANCOVA processing 0 ")
					if(tclvalue(permute.val)=="1" & ncpu>1) cl = makeCluster(ncpu, type = "SOCK")

					if(check) tkmessageBox(title="Analysis of Covariance", message = sprintf("Instead of nested ANCOVA, ANCOVA is applied\nbecause %d samples have only one replicate.",check), icon = "info", type = "ok")

					for(fa in 1:length(Factor_label)){
						fac = Factor_label[fa]
						if(ctrl_label!=""){
							group_idx = which(Factor %in% c(fac, ctrl_label))
							Factor_temp = factor(Factor[group_idx])
						} else {
							group_idx = 1:nrow(Factor)
							Factor_temp = as.numeric(Factor[,fa])
						}
						Rep_temp = factor(Rep[group_idx])
						Ind_temp <- factor(Ind[group_idx],level=unique(Ind[group_idx]))

						options(contrasts = c("contr.helmert", "contr.poly"))
						anova <- apply(peakabun[group_idx,],2,function(x) {
							nLevel <- length(unique(Factor_temp[!is.na(x)]))
							if(any(nLevel<2)) NA else aov(as.formula(ifelse(check,"x ~ Factor_temp","x ~ Factor_temp + Factor_temp:Rep_temp")))
						})
						pvalue <- do.call("rbind",sapply(anova,function(x) if(all(is.na(x))) NA else if(summary(x)[[1]]["Residuals","Sum Sq"]<2.2e-16) NA else Anova(x,type="III",singular.ok=T)[,"Pr(>F)"],simplify=F))
						colnames(pvalue) <- gsub(" ","",gsub("_temp","",gsub("Factor_temp", sprintf("%s%s%s",fac,ifelse(ctrl_label!="","_",""),ctrl_label), c("(Intercept)",rownames(summary(anova[[which(!apply(pvalue,1,function(x) all(is.na(x))))[1]]])[[1]])))))
						options(contrasts = c("contr.SAS", "contr.poly"))
						if(ctrl_label!="") contrasts(Factor_temp) = contrasts(Factor_temp)[c(which(levels(Factor_temp)!=ctrl_label), which(levels(Factor_temp)==ctrl_label)),]
						beta_temp <- apply(peakabun[group_idx,],2,function(x) {
							nLevel <- length(unique(Factor_temp[!is.na(x)]))
							if(any(nLevel<2)) NA else lm(as.formula(ifelse(check,"x ~ Factor_temp","x ~ Factor_temp + Factor_temp:Rep_temp")))$coefficients
						})
						beta_temp = switch(class(beta_temp), "matrix" = t(beta_temp), "list"=do.call(rbind, beta_temp))

						
						
						
						pvalue = pvalue[,-ncol(pvalue)]
						write.table(cbind(peak_index, pvalue), file = sprintf("%s/p-value_%s%s%s_details.csv",outpath,fac,ifelse(ctrl_label!=""," vs. ",""),ctrl_label), col.names = T, row.names = F, quote = 1, sep = ",")
						
						pvalue = as.data.frame(pvalue[,2])
						beta_temp = as.data.frame(beta_temp[,2])
						colnames(pvalue) = sprintf("%s%s%s",fac,ifelse(ctrl_label!="","_",""),ctrl_label)
						pfdr = apply(pvalue, 2, function(x) p.adjust(x, method = "BH"))
						colnames(pvalue) = paste(colnames(pvalue), "(pval)", sep = "")
						colnames(pfdr) = paste(colnames(pfdr), "(pfdr)", sep = "")
						colnames(beta_temp) = sprintf("%s%s%s(coef)",fac,ifelse(ctrl_label!="","_",""),ctrl_label)
						
						write.table(cbind(peak_index, pfdr), file = sprintf("%s/pFDR_%s%s%s.csv",outpath,fac,ifelse(ctrl_label!=""," vs. ",""),ctrl_label), col.names = T, row.names = F, quote = 1, sep = ",")
						write.table(cbind(peak_index, beta_temp), file = sprintf("%s/coef_%s%s%s.csv",outpath,fac,ifelse(ctrl_label!=""," vs. ",""),ctrl_label), col.names = T, row.names = F, quote = 1, sep = ",")
						
						
						
						if(tclvalue(permute.val)=="1"){
							if(sum(pvalue[,1]<as.numeric(tclvalue(permutecut.val)),na.rm=T)>0){
								peak_sig = which(pvalue[,1]<as.numeric(tclvalue(permutecut.val)))
								pvalue_sig = pvalue[peak_sig, 1]
								peakabun_sig = as.data.frame(peakabun[group_idx,peak_sig])
								if(ncpu==1){
									pvalue_permute = sapply(1:as.numeric(tclvalue(permute.times)), function(perm){
										Factor_temp_perm = sample(Factor_temp)
										
										
											
										
											
											
											
										
										options(contrasts = c("contr.helmert", "contr.poly"))

										
										
										anova_temp <- apply(peakabun_sig,2,function(x) {
											nLevel <- length(unique(Factor_temp_perm[!is.na(x)]))
											if(any(nLevel<2)) NA else aov(as.formula(ifelse(check,"x ~ Factor_temp_perm","x ~ Factor_temp_perm + Factor_temp_perm:Rep_temp")))
										})
										pvalue_temp <- do.call("cbind",sapply(anova_temp,function(x) if(all(is.na(x))) NA else if(summary(x)[[1]]["Residuals","Sum Sq"]<2.2e-16) NA else Anova(x,type="III",singular.ok=T)[,"Pr(>F)"],simplify=F))
										info <- sprintf("%d%%", round(100*fa/length(Factor_label)*(perm/as.numeric(tclvalue(permute.times)))))
										setTkProgressBar(pb, value = round(100*fa/length(Factor_label)*(perm/as.numeric(tclvalue(permute.times)))), sprintf("ANCOVA-Please wait for ANCOVA processing (%s done)", info), paste(fac, "-", ctrl_label, info, sep = " "))

										return(pvalue_temp[2,]<pvalue_sig)
									})
								}else if(ncpu>1){
									
									clusterExport(cl, c("check","Factor_temp", "Rep_temp", "peakabun_sig", "pvalue_sig", "Anova", "osIsWin"), envir=environment())
									pvalue_permute = parSapply(cl, 1:as.numeric(tclvalue(permute.times)), function(perm){

										pb <- tkProgressBar(sprintf("ANCOVA-Parallel task (permutation %d)",perm), "0% done", 0, 100, 0, width = 500)
											Factor_temp_perm = sample(Factor_temp)
											
											
												
											
												
												
												
											
											options(contrasts = c("contr.helmert", "contr.poly"))
											
											
											anova_temp <- apply(peakabun_sig,2,function(x) {
												nLevel <- length(unique(Factor_temp_perm[!is.na(x)]))
												if(any(nLevel<2)) NA else aov(as.formula(ifelse(check,"x ~ Factor_temp_perm","x ~ Factor_temp_perm + Factor_temp_perm:Rep_temp")))
											})
										if(osIsWin) for(i in 1:100) setTkProgressBar(pb, value = i, sprintf("ANCOVA-Parallel task (permutation %d)", perm), "")

										pvalue_temp <- do.call("cbind",sapply(anova_temp,function(x) if(all(is.na(x))) NA else if(summary(x)[[1]]["Residuals","Sum Sq"]<2.2e-16) NA else Anova(x,type="III",singular.ok=T)[,"Pr(>F)"],simplify=F))
										return(pvalue_temp[2,]<pvalue_sig)
									})
								}
								empvalue = as.data.frame((rowSums(pvalue_permute)+1)/(as.numeric(tclvalue(permute.times))+1))
								empvalue = cbind(empvalue, p.adjust(c(empvalue[,1], rep(1, nrow(peak_index)-nrow(empvalue))), method = "BH")[1:nrow(empvalue)])
								colnames(empvalue) = c(paste(fac, "_", ctrl_label, "(epv)", sep = ""), paste(fac, "_", ctrl_label, "(epFDR)", sep = ""))
								write.table(cbind(peak_index[peak_sig,], empvalue), file = sprintf("%s/Empirical p-value_%s%s%s.csv",outpath,fac,ifelse(ctrl_label!=""," vs. ",""),ctrl_label), col.names = T, row.names = F, quote = 1, sep = ",")
							}else{
								tkmessageBox(title = "Warning", message = "No significant peaks", icon = "warning", type = "ok")
							}
						}
						

						info <- sprintf("%d%%", round(100*fa/length(Factor_label)))
						setTkProgressBar(pb, value = round(100*fa/length(Factor_label)), sprintf("ANCOVA-Please wait for ANCOVA processing (%s done)", info), paste(fac, "-", ctrl_label, info, sep = " "))
						if(round(100*fa/length(Factor_label))<100){
							cat(round(100*fa/length(Factor_label)), " ", sep = "")
						}else{
							cat(round(100*fa/length(Factor_label)), " \n", sep = "")
						}
						Sys.sleep(0.1)
					}
					if(tclvalue(permute.val)=="1" & ncpu>1) stopCluster(cl)
					options(op)
					setTkProgressBar(pb, value = 100, "ANCOVA-Please wait for ANCOVA processing (100% done)", "Finished 100%")
					Sys.sleep(1)
					close(pb)
					tkconfigure(tt,cursor="arrow")
					
					toppeak = tktoplevel(); if(isIcon) tk2ico.set(toppeak,icon)
					tkwm.title(toppeak, "Top peaks")
					fr_top = tkframe(toppeak)
					tkgrid(fr_top)
					top.lab <- tklabel(fr_top, text = "    How many top peaks do you want to paint red?  ")
					top.num = tclVar("50")
					top.entry <- tkentry(fr_top, width = 6, textvariable = top.num, bg = "white", validatecommand = "string is integer %P", validate = "all")
					
					tkgrid(top.lab, top.entry, tklabel(fr_top, text = paste("(1 ~ ", ncol(peakabun), ")     ", sep = "")))
					tkgrid(tklabel(toppeak, text = "", height = 0, font = fontIntro_para))
					fr_toppeak = tkframe(toppeak)
					tkgrid(fr_toppeak)
					plot.but <- tkbutton(fr_toppeak, text = "  Plot  ", command = function(...){
						tkconfigure(tt,cursor="watch")
						top.num = as.numeric(tclvalue(top.num))
						if(!top.num%in%seq.int(1, ncol(peakabun))) top.num = 50
						for(fa in 1:length(Factor_label)){
							fac = Factor_label[fa]
							pfdr = read.table(sprintf("%s/pFDR_%s%s%s.csv",outpath,fac,ifelse(ctrl_label!=""," vs. ",""),ctrl_label), sep = ",", header = T, quote = "\"")
							pvalue = read.table(sprintf("%s/p-value_%s%s%s_details.csv",outpath,fac,ifelse(ctrl_label!=""," vs. ",""),ctrl_label), sep = ",", header = T, quote = "\"")
							beta_temp = read.table(sprintf("%s/coef_%s%s%s.csv",outpath,fac,ifelse(ctrl_label!=""," vs. ",""),ctrl_label), sep = ",", header = T, quote = "\"")
							

							png(filename = sprintf("%s/Volcanoplot_%s%s%s_pFDR_%d.png",outpath,fac,ifelse(ctrl_label!=""," vs. ",""),ctrl_label,top.num), width = 1440, height = 800)
							i=4
							par(mar = c(5,4,1,3))
							pch = rep(21, nrow(pfdr))
							pch[which(-log10(pfdr[,i])>=-log10(0.05))] = 24
							bg = rep("white", nrow(pfdr))
							bg[which(-log10(pfdr[,i])>=-log10(0.05))] = "green3"
							bg[order(-log10(pfdr[,i]), decreasing = T)[1:top.num]] = "red"
							plot(beta_temp[,4], -log10(pfdr[,i]), xlab = "", ylab = "", cex = 1, pch = pch, bg = bg)
							text(beta_temp[order(-log10(pfdr[,i]), decreasing = T)[1:top.num],4], -log10(pfdr[order(-log10(pfdr[,i]), decreasing = T)[1:top.num],i]), label = peak_index[order(-log10(pfdr[,i]), decreasing = T)[1:top.num],1], cex = 1, pos = 1, adj = 0.2)

							abline(h = -log10(0.05), col = "green3")
							par(xpd = NA)
							text(par("usr")[1], -log10(0.05), col = "green3", label = substitute(paste(plain(-log[10]), "(0.05) (", test, " significant peaks)"), list(test = sum(pfdr[,i]<0.05))), font = 4, cex = 1.5, adj = c(0,-0.1))
							par(xpd = F)
							mtext("Fold Change (coefficient)", line = 2, side = 1, col="black", padj = 1, cex = 1.5)
							mtext(expression(-log[10](pFDR)), line = 3.5, side = 2, col="black", padj = 1, cex = 1.5)
							legend(x = "topleft", legend = c("p-FDR<0.05", paste("top", top.num, sep = "")), pch = 24, bty = "o", box.col = "black", cex = 1, pt.bg = c("green3", "red"), bg = rgb(1, 1, 1, alpha = 0.7))
							dev.off()

							png(filename = sprintf("%s/Volcanoplot_%s%s%s_pv_%d.png",outpath,fac,ifelse(ctrl_label!=""," vs. ",""),ctrl_label,top.num), width = 1440, height = 800)
							i=5
							par(mar = c(5,4,1,3))
							pch = rep(21, nrow(pvalue))
							pch[which(-log10(pvalue[,i])>=-log10(0.05))] = 24
							bg = rep("white", nrow(pvalue))
							bg[which(-log10(pvalue[,i])>=-log10(0.05))] = "green3"
							bg[order(-log10(pvalue[,i]), decreasing = T)[1:top.num]] = "red"
							plot(beta_temp[,4], -log10(pvalue[,i]), cex = 1, pch = pch, bg = bg, xlab = "", ylab = "")
							text(beta_temp[order(-log10(pvalue[,i]), decreasing = T)[1:top.num],4], -log10(pvalue[order(-log10(pvalue[,i]), decreasing = T)[1:top.num],i]), label = peak_index[order(-log10(pvalue[,i]), decreasing = T)[1:top.num],1], cex = 1, pos = 1, adj = 0.2)
							abline(h = -log10(0.05), col = "green3")
							par(xpd = NA)
							text(par("usr")[1], -log10(0.05), col = "green3", label = substitute(paste(plain(-log[10]), "(0.05) (", test, " significant peaks)"), list(test = sum(pvalue[,i]<0.05))), font = 4, cex = 1.5, adj = c(0,-0.1))
							par(xpd = F)
							mtext("Fold Change (coefficient)", line = 2, side = 1, col="black", padj = 1, cex = 1.5)
							mtext(expression(-log[10](pv)), line = 3.5, side = 2, col="black", padj = 1, cex = 1.5)
							legend(x = "topleft", legend = c("p-value<0.05", paste("top", top.num, sep = "")), pch = 24, bty = "o", box.col = "black", cex = 1, pt.bg = c("green3", "red"), bg = rgb(1, 1, 1, alpha = 0.7))
							dev.off()
							

							
							if(tclvalue(permute.val)=="1" & file.exists(sprintf("%s/Empirical p-value_%s%s%s.csv",outpath,fac,ifelse(ctrl_label!=""," vs. ",""),ctrl_label))){
								beta_temp = beta_temp[which(pvalue[,5]<as.numeric(tclvalue(permutecut.val))),]
								peak_index_temp = peak_index[which(pvalue[,5]<as.numeric(tclvalue(permutecut.val))),]
								pvalue = read.table(sprintf("%s/Empirical p-value_%s%s%s.csv",outpath,fac,ifelse(ctrl_label!=""," vs. ",""),ctrl_label), sep = ",", header = T, quote = "\"")
								pfdr = pvalue
								


								png(filename = sprintf("%s/Volcanoplot_%s%s%s_epFDR_%d.png",outpath,fac,ifelse(ctrl_label!=""," vs. ",""),ctrl_label,top.num), width = 1440, height = 800)
								i=5
								par(mar = c(5,4,1,3))
								pch = rep(21, nrow(pfdr))
								pch[which(-log10(pfdr[,i])>=-log10(0.05))] = 24
								bg = rep("white", nrow(pfdr))
								bg[which(-log10(pfdr[,i])>=-log10(0.05))] = "green3"
								bg[order(-log10(pfdr[,i]), decreasing = T)[1:min(top.num, nrow(pfdr))]] = "red"
								plot(beta_temp[,4], -log10(pfdr[,i]), xlab = "", ylab = "", cex = 1, pch = pch, bg = bg)
								text(beta_temp[order(-log10(pfdr[,i]), decreasing = T)[1:min(top.num, nrow(pfdr))],4], -log10(pfdr[order(-log10(pfdr[,i]), decreasing = T)[1:min(top.num, nrow(pfdr))],i]), label = peak_index_temp[order(-log10(pfdr[,i]), decreasing = T)[1:min(top.num, nrow(pfdr))],1], cex = 1, pos = 1, adj = 0.2)

								abline(h = -log10(0.05), col = "green3")
								par(xpd = NA)
								text(par("usr")[1], -log10(0.05), col = "green3", label = substitute(paste(plain(-log[10]), "(0.05) (", test, " significant peaks)"), list(test = sum(pfdr[,i]<0.05))), font = 4, cex = 1.5, adj = c(0,-0.1))
								par(xpd = F)
								mtext("Fold Change (coefficient)", line = 2, side = 1, col="black", padj = 1, cex = 1.5)
								mtext(expression(-log[10](epFDR)), line = 3.5, side = 2, col="black", padj = 1, cex = 1.5)
								legend(x = "topleft", legend = c("Empirical p-FDR<0.05", paste("top", min(top.num, nrow(pfdr)), sep = "")), pch = 24, bty = "o", box.col = "black", cex = 1, pt.bg = c("green3", "red"), bg = rgb(1, 1, 1, alpha = 0.7))
								dev.off()
								
								png(filename = sprintf("%s/Volcanoplot_%s%s%s_epv_%d.png",outpath,fac,ifelse(ctrl_label!=""," vs. ",""),ctrl_label,top.num), width = 1440, height = 800)
								i=4
								par(mar = c(5,4,1,3))
								pch = rep(21, nrow(pvalue))
								pch[which(-log10(pvalue[,i])>=-log10(0.05))] = 24
								bg = rep("white", nrow(pvalue))
								bg[which(-log10(pvalue[,i])>=-log10(0.05))] = "green3"
								bg[order(-log10(pvalue[,i]), decreasing = T)[1:min(top.num, nrow(pvalue))]] = "red"
								plot(beta_temp[,4], -log10(pvalue[,i]), cex = 1, pch = pch, bg = bg, xlab = "", ylab = "")
								text(beta_temp[order(-log10(pvalue[,i]), decreasing = T)[1:min(top.num, nrow(pvalue))],4], -log10(pvalue[order(-log10(pvalue[,i]), decreasing = T)[1:min(top.num, nrow(pvalue))],i]), label = peak_index_temp[order(-log10(pvalue[,i]), decreasing = T)[1:min(top.num, nrow(pvalue))],1], cex = 1, pos = 1, adj = 0.2)
								abline(h = -log10(0.05), col = "green3")
								par(xpd = NA)
								text(par("usr")[1], -log10(0.05), col = "green3", label = substitute(paste(plain(-log[10]), "(0.05) (", test, " significant peaks)"), list(test = sum(pvalue[,i]<0.05))), font = 4, cex = 1.5, adj = c(0,-0.1))
								par(xpd = F)
								mtext("Fold Change (coefficient)", line = 2, side = 1, col="black", padj = 1, cex = 1.5)
								mtext(expression(-log[10](epv)), line = 3.5, side = 2, col="black", padj = 1, cex = 1.5)
								legend(x = "topleft", legend = c("Empirical p-value<0.05", paste("top", min(top.num, nrow(pvalue)), sep = "")), pch = 24, bty = "o", box.col = "black", cex = 1, pt.bg = c("green3", "red"), bg = rgb(1, 1, 1, alpha = 0.7))
								dev.off()

							}
						}
						tkmessageBox(title = "Analysis of Covariance", message = "Graphs are saved.", icon = "info", type = "ok")
						cat("ANCOVA - Graphs are saved.\n", sep = "")
						tkconfigure(tt,cursor="arrow")
					})
					close.but <- tkbutton(fr_toppeak, text = "  Close  ", command = function(...){
						tkdestroy(toppeak)
						tkmessageBox(title="Analysis of Covariance", message = "ANCOVA is done.", icon = "info", type = "ok")
						cat("ANCOVA - ANCOVA is done.\n", sep = "")
					})
					
					tkgrid(plot.but, tklabel(fr_toppeak, text = "                              "), close.but)
					tkgrid(tklabel(toppeak, text = "", height = 0, font = fontIntro_para))
				}else{
					tkmessageBox(title = "Error", message = "Sample ID not match.\nPlease input the adequate factors file.", icon = "error", type = "ok")
					cat("ANCOVA(Error) - Sample ID not match. Please input the adequate factors file.\n")
					tkfocus(dlg)
				}
			}else if(!(tclvalue(textbatchinput)=="") & (tclvalue(textcova)=="")){  
				batch_effect = read.table(tclvalue(textbatchinput), header = T, sep = ifelse(grepl(".txt", tclvalue(textbatchinput)), "\t", ","), as.is = T, na.string = c("NA", as.character(tclvalue(text_missing))), quote = "\"")
				
				peak_index = peakabun[,c(1:3)]
				peakabun = peakabun[,-c(1:3)]
				colnames(peakabun) = gsub("^\\d{4,8}_(.*?)|(^[A-Z]\\d{4,8})_(.*?)", "\\1", colnames(peakabun))
				Factor[,1] = gsub("^\\d{4,8}_(.*?)|(^[A-Z]\\d{4,8})_(.*?)", "\\1", Factor[,1])
				batch_effect[,1] = gsub("^\\d{4,8}_(.*?)|(^[A-Z]\\d{4,8})_(.*?)", "\\1", batch_effect[,1])
				
				if(all(!is.na(match(gsub("_[0-9]", "", colnames(peakabun)), Factor[,1])))& all(!is.na(match(colnames(peakabun), batch_effect[,1]))))
				{
					tkdestroy(dlg)
					tkconfigure(tt,cursor="watch")
					
					Rep = gsub("(.*?)_(.*?)", "\\2", colnames(peakabun))
					Ind <- gsub("_.*", "", colnames(peakabun))

					if(ctrl_label!=""){
						Factor = Factor[match(gsub("_[0-9]", "", colnames(peakabun)), Factor[,1]),2]
						Factor_label = levels(factor(Factor))
						Factor_label = Factor_label[-which(Factor_label==ctrl_label)]
					} else {
						Factor <- Factor[match(gsub("_[0-9]", "", colnames(peakabun)), Factor[,1]),]
						Factor <- subset(Factor,select=-1) 
						Factor_label <- colnames(Factor)
					}

					batch_effect = if(ncol(batch_effect)>2){
						as.data.frame(batch_effect[match(colnames(peakabun), batch_effect[,1]),-1])
					}else{
						if(grepl("latent_group_", tclvalue(textbatchinput))){
							as.data.frame(matrix(factor(batch_effect[match(colnames(peakabun), batch_effect[,1]),-1]), ncol = 1, dimnames = list(NULL, colnames(batch_effect)[2])))
						} else {
							as.data.frame(matrix(batch_effect[match(colnames(peakabun), batch_effect[,1]),-1], ncol = 1, dimnames = list(NULL, colnames(batch_effect)[2])))
						}
					}
					

					peakabun = t(peakabun)
					if(ctrl_label!=""){
						Factor_label = levels(factor(Factor))
						Factor_label = Factor_label[-which(Factor_label==ctrl_label)]
					} else {
						Factor_label <- colnames(Factor)
					}
					
					
					fontHeading <- tkfont.create(family="times",size=15,slant="italic")
					
					vif_check <- tktoplevel(); if(isIcon) tk2ico.set(vif_check,icon)
					tkwm.title(vif_check, "Analysis of Covariance - Variance inflation Factor (VIF)")
					tkpack(fr_head <- tkframe(vif_check), side = "top")
					tkpack(tklabel(fr_head, text = paste("Collinearity Check - VIF", sep = ""), font = fontHeading), side = "top")
					tkpack(fr_vif <- tkframe(vif_check), side = "top")
					
					tclarray <- tclArray()
					VIF_table = c(ifelse(ctrl_label!="","Group_Comparison",""), "Variable", colnames(batch_effect))
					
					for(lab in Factor_label){
						if(ctrl_label!="") case_idx = which(Factor %in% c(lab, ctrl_label))
						else case_idx = 1:nrow(Factor)
						temp_factor = c(ifelse(ctrl_label!="",paste(lab, "vs.", ctrl_label, sep = ""),lab), ifelse(ctrl_label!="",paste(lab, "-", ctrl_label, "(var)", sep = ""),paste(lab, "(var)", sep = "")), 
						sapply(1:ncol(batch_effect), function(x){
							if(ctrl_label!="") temp = vif(lm(peakabun[case_idx,1]~factor(Factor[case_idx]) + batch_effect[case_idx,x]))
							else temp = vif(lm(peakabun[case_idx,1]~Factor[case_idx,lab] + batch_effect[case_idx,x]))
							ifelse(class(temp)=="matrix", round(temp[2,3], 1), round(temp[2], 1))
						}))

						VIF_table = cbind(VIF_table, temp_factor)
					}
					
					for (i in (1:nrow(VIF_table)))
						for (j in (1:ncol(VIF_table)))
							tclarray[[i-1,j-1]] <- VIF_table[i,j]
					scr <- tkscrollbar(fr_vif, repeatinterval = 35, orient="horizontal", command = function(...) tkxview(table1,...))
					table1 <- tkwidget(fr_vif,"table",variable=tclarray,rows=nrow(VIF_table),cols=ncol(VIF_table),titlerows=2,titlecols = 1, selectmode="extended",colwidth=20,background="white", xscrollcommand = function(...) tkset(scr, ...))
					tkpack(table1)
					tkpack(scr, fill = "x")
					tkconfigure(table1,selectmode="extended",rowseparator="\"\n\"",colseparator="\"\t\"")
					tkconfigure(table1,resizeborders="none", state = "disable")    
					tkpack(tklabel(fr_vif, text = "", height = 0, font = fontIntro_para), side = "bottom")

					tkpack(tklabel(vif_check, text = "", height = 0, font = fontIntro_para), side = "bottom")
					tkpack(fr_next <- tkframe(vif_check), side = "bottom", fill = "x")
					vif_threshold <- tclVar("10")
					tkpack(tklabel(fr_next,text="       The VIF upper bound:    "), side = "left")
					tkpack(tkentry(fr_next, width="7", textvariable = vif_threshold, bg = "white"), side = "left")
					tkpack(tklabel(fr_next, text = "  ", width = 3), side = "right")
					

					
					tkpack(tkbutton(fr_next, text = "  Next  ", command = function(...){
						tkdestroy(vif_check)
						
						tkconfigure(tt,cursor="watch")
						
						vif_threshold = as.numeric(tclvalue(vif_threshold))
						write.table(VIF_table, file = paste(outpath, "/VIF_threshold_", vif_threshold, ".csv", sep = ""), sep = ",", col.names = F, row.names = F, quote = F)
						VIF_table = matrix(as.numeric(VIF_table[-c(1,2),-1]), ncol = ncol(VIF_table)-1)
						VIF_table1 = matrix(1, nrow = nrow(VIF_table), ncol = ncol(VIF_table))
						VIF_table1[VIF_table > vif_threshold] = 0
						
						VIF_table = VIF_table1

						if(check) tkmessageBox(title="Analysis of Covariance", message = sprintf("Instead of nested ANCOVA, ANCOVA is applied\nbecause %d samples have only one replicate.",check), icon = "info", type = "ok")
						
						op <- options(contrasts = c("contr.helmert", "contr.poly"))
						pb <- tkProgressBar("ANCOVA - Please wait for ANCOVA processing", "0% done", 0, 100, 0, width = 500)
						cat("ANCOVA - Please wait for ANCOVA processing 0 ")
						if(tclvalue(permute.val)=="1" & ncpu>1) cl = makeCluster(ncpu, type = "SOCK")
						for(fa in 1:length(Factor_label)){
							fac = Factor_label[fa]
							if(ctrl_label!=""){
								group_idx = which(Factor %in% c(fac, ctrl_label))
								Factor_temp = factor(Factor[group_idx])
							} else {
								group_idx <- 1:nrow(Factor)
								Factor_temp <- Factor[,fa]
							}
							Rep_temp = factor(Rep[group_idx])
							Ind_temp = factor(Ind[group_idx],level=unique(Ind[group_idx]))

							be_choose = as.logical(as.numeric(VIF_table[,1]))
							temp_all = as.data.frame(batch_effect[group_idx, which(be_choose)])
							colnames(temp_all) = colnames(batch_effect)[which(be_choose)]


							options(contrasts = c("contr.helmert", "contr.poly"))
							anova <- apply(peakabun[group_idx,],2,function(x) {
								nLevel <- apply(cbind(Factor_temp,temp_all)[!is.na(x),],2,function(y) length(unique(y)))
								if(any(nLevel<2)) NA else aov(as.formula(paste(ifelse(check,"x ~ Factor_temp + ","x ~ Factor_temp + Factor_temp:Rep_temp + "), paste(colnames(temp_all), collapse = "+"))), data = temp_all)
							})
							pvalue <- do.call("rbind",sapply(anova,function(x) if(all(is.na(x))) NA else if(summary(x)[[1]]["Residuals","Sum Sq"]<2.2e-16) NA else Anova(x,type="III",singular.ok=T)[,"Pr(>F)"],simplify=F))
							colnames(pvalue) <- gsub(" ","",gsub("_temp","",gsub("Factor_temp", sprintf("%s%s%s",fac,ifelse(ctrl_label!="","_",""),ctrl_label), c("(Intercept)",rownames(summary(anova[[which(!apply(pvalue,1,function(x) all(is.na(x))))[1]]])[[1]])))))

							options(contrasts = c("contr.SAS", "contr.poly"))
							if(ctrl_label!="") contrasts(Factor_temp) = contrasts(Factor_temp)[c(which(levels(Factor_temp)!=ctrl_label), which(levels(Factor_temp)==ctrl_label)),]
							beta_temp <- apply(peakabun[group_idx,],2,function(x) {
								nLevel <- apply(cbind(Factor_temp,temp_all)[!is.na(x),],2,function(y) length(unique(y)))
								if(any(nLevel<2)) NA else lm(as.formula(paste(ifelse(check,"x ~ Factor_temp + ","x ~ Factor_temp + Factor_temp:Rep_temp + "), paste(colnames(temp_all), collapse = "+"))), data = temp_all)$coefficients
							})
							beta_temp = switch(class(beta_temp), "matrix" = t(beta_temp), "list"=do.call(rbind, beta_temp))
							resi_temp = apply(peakabun[group_idx,], 2, function(x){
								nLevel <- apply(cbind(Factor_temp,temp_all)[!is.na(x),],2,function(y) length(unique(y)))
								temp <- rep(NA, length(x))
								if(all(nLevel>=2)) temp[!is.na(x)&rowSums(is.na(temp_all))==0] <- lm(as.formula(paste(ifelse(check,"x ~ ","x ~ Factor_temp:Rep_temp + "), paste(colnames(temp_all), collapse = "+"))), data = temp_all)$residuals
								return(temp)
							})
							resi_temp = switch(class(resi_temp), "matrix" = t(resi_temp), "list"=do.call(rbind, resi_temp))

							
							colnames(pvalue) <- gsub(" ","",gsub("_temp","",gsub("Factor_temp", sprintf("%s%s%s",fac,ifelse(ctrl_label!="","_",""),ctrl_label), c("(Intercept)",rownames(summary(anova[[1]])[[1]])))))
							
							
							

							pvalue = pvalue[,-ncol(pvalue)]
							write.table(cbind(peak_index, pvalue), file = sprintf("%s/p-value_%s%s%s_details.csv",outpath,fac,ifelse(ctrl_label!=""," vs. ",""),ctrl_label), col.names = T, row.names = F, quote = 1, sep = ",")
							pvalue = as.data.frame(pvalue[,2])
							beta_temp = as.data.frame(beta_temp[,2])
							resi_temp = as.data.frame(resi_temp); colnames(resi_temp) <- rownames(peakabun[group_idx,]) 
							colnames(pvalue) = sprintf("%s%s%s",fac,ifelse(ctrl_label!="","_",""),ctrl_label)
							pfdr = apply(pvalue, 2, function(x) p.adjust(x, method = "BH"))
							colnames(pvalue) = paste(colnames(pvalue), "(pval)", sep = "")
							colnames(pfdr) = paste(colnames(pfdr), "(pfdr)", sep = "")
							colnames(beta_temp) = sprintf("%s%s%s(coef)",fac,ifelse(ctrl_label!="","_",""),ctrl_label)
							
							write.table(cbind(peak_index, pfdr), file = sprintf("%s/pFDR_%s%s%s.csv",outpath,fac,ifelse(ctrl_label!=""," vs. ",""),ctrl_label), col.names = T, row.names = F, quote = 1, sep = ",")
							write.table(cbind(peak_index, beta_temp), file = sprintf("%s/coef_%s%s%s.csv",outpath,fac,ifelse(ctrl_label!=""," vs. ",""),ctrl_label), col.names = T, row.names = F, quote = 1, sep = ",")
							write.table(cbind(peak_index, resi_temp), file = sprintf("%s/resi_%s%s%s_forPLS.csv",outpath,fac,ifelse(ctrl_label!=""," vs. ",""),ctrl_label), col.names = T, row.names = F, quote = 1, sep = ",")
							if(tclvalue(permute.val)!="1") write.table(cbind(peak_index, resi_temp)[pfdr<.05,], file = sprintf("%s/resi_pFDRsig_%s%s%s_forPLS.csv",outpath,fac,ifelse(ctrl_label!=""," vs. ",""),ctrl_label), col.names = T, row.names = F, quote = 1, sep = ",")
							
							if(tclvalue(permute.val)=="1"){
								if(sum(pvalue[,1]<as.numeric(tclvalue(permutecut.val)),na.rm=T)>0){
									peak_sig = which(pvalue[,1]<as.numeric(tclvalue(permutecut.val)))
									pvalue_sig = pvalue[peak_sig, 1]
									peakabun_sig = as.data.frame(peakabun[group_idx,peak_sig])
									if(ncpu==1){
										pvalue_permute = sapply(1:as.numeric(tclvalue(permute.times)), function(perm){
											Factor_temp_perm = sample(Factor_temp)
											
											
												
											
												
												
												
											
											options(contrasts = c("contr.helmert", "contr.poly"))

											
											
											
											
											anova_temp <- apply(peakabun_sig,2,function(x) {
												nLevel <- length(unique(Factor_temp_perm[!is.na(x)]))
												if(any(nLevel<2)) NA else aov(as.formula(paste(ifelse(check,"x ~ Factor_temp_perm","x ~ Factor_temp_perm + Factor_temp_perm:Rep_temp + "),paste(colnames(temp_all), collapse = "+"))),data = temp_all)
											})
											pvalue_temp <- do.call("cbind",sapply(anova_temp,function(x) if(all(is.na(x))) NA else if(summary(x)[[1]]["Residuals","Sum Sq"]<2.2e-16) NA else Anova(x,type="III",singular.ok=T)[,"Pr(>F)"],simplify=F))
											info <- sprintf("%d%%", round(100*fa/length(Factor_label)*(perm/as.numeric(tclvalue(permute.times)))))
											setTkProgressBar(pb, value = round(100*fa/length(Factor_label)*(perm/as.numeric(tclvalue(permute.times)))), sprintf("ANCOVA-Please wait for ANCOVA processing (%s done)", info), paste(fac, "-", ctrl_label, info, sep = " "))

											return(pvalue_temp[2,]<pvalue_sig)
										})
									}else if(ncpu>1){
										clusterExport(cl, c("check","Factor_temp", "Rep_temp", "peakabun_sig", "pvalue_sig", "temp_all", "Anova", "osIsWin"), envir=environment())
										
										pvalue_permute = parSapply(cl, 1:as.numeric(tclvalue(permute.times)), function(perm){

											pb <- tkProgressBar(sprintf("ANCOVA-Parallel task (permutation %d)",perm), "0% done", 0, 100, 0, width = 500)
												Factor_temp_perm = sample(Factor_temp)
												
												
													
												
													
													
													
												
												options(contrasts = c("contr.helmert", "contr.poly"))

												
												
												
												
												anova_temp <- apply(peakabun_sig,2,function(x) {
													nLevel <- length(unique(Factor_temp_perm[!is.na(x)]))
													if(any(nLevel<2)) NA else aov(as.formula(paste(ifelse(check,"x ~ Factor_temp_perm","x ~ Factor_temp_perm + Factor_temp_perm:Rep_temp + "),paste(colnames(temp_all), collapse = "+"))),data = temp_all)
												})
											if(osIsWin) for(i in 1:100) setTkProgressBar(pb, value = i, sprintf("ANCOVA-Parallel task (permutation %d)", perm), "")

											pvalue_temp <- do.call("cbind",sapply(anova_temp,function(x) if(all(is.na(x))) NA else if(summary(x)[[1]]["Residuals","Sum Sq"]<2.2e-16) NA else Anova(x,type="III",singular.ok=T)[,"Pr(>F)"],simplify=F))
											return(pvalue_temp[2,]<pvalue_sig)
										})
									}
									empvalue = as.data.frame((rowSums(pvalue_permute)+1)/(as.numeric(tclvalue(permute.times))+1))
									empvalue = cbind(empvalue, p.adjust(c(empvalue[,1], rep(max(c(empvalue[,1], as.numeric(tclvalue(permutecut.val)))), nrow(peak_index)-nrow(empvalue))), method = "BH")[1:nrow(empvalue)])
									colnames(empvalue) = c(sprintf("%s%s%s(epv)",fac,ifelse(ctrl_label!="","_",""),ctrl_label), sprintf("%s%s%s(epFDR)",fac,ifelse(ctrl_label!="","_",""),ctrl_label))
									write.table(cbind(peak_index[peak_sig,], empvalue), file = sprintf("%s/Empirical p-value_%s%s%s.csv",outpath,fac,ifelse(ctrl_label!=""," vs. ",""),ctrl_label), col.names = T, row.names = F, quote = 1, sep = ",")
									write.table(cbind(peak_index, resi_temp)[peak_sig,][empvalue[,2]<.05,], file = sprintf("%s/resi_EmpiricalpFDRsig_%s%s%s_forPLS.csv",outpath,fac,ifelse(ctrl_label!=""," vs. ",""),ctrl_label), col.names = T, row.names = F, quote = 1, sep = ",")
								}else{
									tkmessageBox(title = "Warning", message = "No significant peaks", icon = "warning", type = "ok")
								}
							}
						
							
							VIF_table = matrix(VIF_table[,-1], nrow = nrow(VIF_table))
							info <- sprintf("%d%%", round(100*fa/length(Factor_label)))
							setTkProgressBar(pb, value = round(100*fa/length(Factor_label)), sprintf("ANCOVA - Please wait for ANCOVA processing (%s done)", info), paste(fac, "-", ctrl_label, info, sep = " "))
							if(round(100*fa/length(Factor_label))<100){
								cat(round(100*fa/length(Factor_label)), " ", sep = "")
							}else{
								cat(round(100*fa/length(Factor_label)), " \n", sep = "")
							}
							Sys.sleep(0.1)
						}
						if(tclvalue(permute.val)=="1" & ncpu>1) stopCluster(cl)
						setTkProgressBar(pb, value = 100, "ANCOVA - Please wait for ANCOVA processing (100% done)", "Finished 100%")
						Sys.sleep(1)
						close(pb)
						options(op)
						tkconfigure(tt,cursor="arrow")
						
						
						toppeak = tktoplevel(); if(isIcon) tk2ico.set(toppeak,icon)
						tkwm.title(toppeak, "Top peaks")
						fr_top = tkframe(toppeak)
						tkgrid(fr_top)
						top.lab <- tklabel(fr_top, text = "    How many top peaks do you want to paint red?  ")
						top.num = tclVar("50")
						top.entry <- tkentry(fr_top, width = 6, textvariable = top.num, bg = "white",validatecommand="string is integer %P", validate = "all")
						
						tkgrid(top.lab, top.entry, tklabel(fr_top, text = paste("(1 ~ ", ncol(peakabun), ")     ", sep = "")))
						tkgrid(tklabel(toppeak, text = "", height = 0, font = fontIntro_para))
						fr_toppeak = tkframe(toppeak)
						tkgrid(fr_toppeak)
						plot.but <- tkbutton(fr_toppeak, text = "  Plot  ", command = function(...){
							tkconfigure(tt,cursor="watch")
							top.num = as.numeric(tclvalue(top.num))
							if(!top.num%in%seq.int(1, ncol(peakabun))) top.num = 50
							for(fa in 1:length(Factor_label)){
								pfdr = read.table(sprintf("%s/pFDR_%s%s%s.csv",outpath,fac,ifelse(ctrl_label!=""," vs. ",""),ctrl_label), sep = ",", header = T, quote = "\"")
								pvalue = read.table(sprintf("%s/p-value_%s%s%s_details.csv",outpath,fac,ifelse(ctrl_label!=""," vs. ",""),ctrl_label), sep = ",", header = T, quote = "\"")
								beta_temp = read.table(sprintf("%s/coef_%s%s%s.csv",outpath,fac,ifelse(ctrl_label!=""," vs. ",""),ctrl_label), sep = ",", header = T, quote = "\"")


								png(filename = sprintf("%s/Volcanoplot_%s%s%s_pFDR_%d.png",outpath,fac,ifelse(ctrl_label!=""," vs. ",""),ctrl_label,top.num), width = 1440, height = 800)
								i=4
								par(mar = c(5,4,1,3))
								pch = rep(21, nrow(pfdr))
								pch[which(-log10(pfdr[,i])>=-log10(0.05))] = 24
								bg = rep("white", nrow(pfdr))
								bg[which(-log10(pfdr[,i])>=-log10(0.05))] = "green3"
								bg[order(-log10(pfdr[,i]), decreasing = T)[1:top.num]] = "red"
								plot(beta_temp[,4], -log10(pfdr[,i]), xlab = "", ylab = "", cex = 1, pch = pch, bg = bg)
								text(beta_temp[order(-log10(pfdr[,i]), decreasing = T)[1:top.num],4], -log10(pfdr[order(-log10(pfdr[,i]), decreasing = T)[1:top.num],i]), label = peak_index[order(-log10(pfdr[,i]), decreasing = T)[1:top.num],1], cex = 1, pos = 1, adj = 0.2)

								abline(h = -log10(0.05), col = "green3")
								par(xpd = NA)
								text(par("usr")[1], -log10(0.05), col = "green3", label = substitute(paste(plain(-log[10]), "(0.05) (", test, " significant peaks)"), list(test = sum(pfdr[,i]<0.05))), font = 4, cex = 1.5, adj = c(0,-0.1))
								par(xpd = F)
								mtext("Fold Change (coefficient)", line = 2, side = 1, col="black", padj = 1, cex = 1.5)
								mtext(expression(-log[10](pFDR)), line = 3.5, side = 2, col="black", padj = 1, cex = 1.5)
								legend(x = "topleft", legend = c("p-FDR<0.05", paste("top", top.num, sep = "")), pch = 24, bty = "o", box.col = "black", cex = 1, pt.bg = c("green3", "red"), bg = rgb(1, 1, 1, alpha = 0.7))
								dev.off()
								
								png(filename = sprintf("%s/Volcanoplot_%s%s%s_pv_%d.png",outpath,fac,ifelse(ctrl_label!=""," vs. ",""),ctrl_label,top.num), width = 1440, height = 800)
								i=5
								par(mar = c(5,4,1,3))
								pch = rep(21, nrow(pvalue))
								pch[which(-log10(pvalue[,i])>=-log10(0.05))] = 24
								bg = rep("white", nrow(pvalue))
								bg[which(-log10(pvalue[,i])>=-log10(0.05))] = "green3"
								bg[order(-log10(pvalue[,i]), decreasing = T)[1:top.num]] = "red"
								plot(beta_temp[,4], -log10(pvalue[,i]), cex = 1, pch = pch, bg = bg, xlab = "", ylab = "")
								text(beta_temp[order(-log10(pvalue[,i]), decreasing = T)[1:top.num],4], -log10(pvalue[order(-log10(pvalue[,i]), decreasing = T)[1:top.num],i]), label = peak_index[order(-log10(pvalue[,i]), decreasing = T)[1:top.num],1], cex = 1, pos = 1, adj = 0.2)
								abline(h = -log10(0.05), col = "green3")
								par(xpd = NA)
								text(par("usr")[1], -log10(0.05), col = "green3", label = substitute(paste(plain(-log[10]), "(0.05) (", test, " significant peaks)"), list(test = sum(pvalue[,i]<0.05))), font = 4, cex = 1.5, adj = c(0,-0.1))
								par(xpd = F)
								mtext("Fold Change (coefficient)", line = 2, side = 1, col="black", padj = 1, cex = 1.5)
								mtext(expression(-log[10](pv)), line = 3.5, side = 2, col="black", padj = 1, cex = 1.5)
								legend(x = "topleft", legend = c("p-value<0.05", paste("top", top.num, sep = "")), pch = 24, bty = "o", box.col = "black", cex = 1, pt.bg = c("green3", "red"), bg = rgb(1, 1, 1, alpha = 0.7))
								dev.off()
								

								
								if(tclvalue(permute.val)=="1" & file.exists(sprintf("%s/Empirical p-value_%s%s%s.csv",outpath,fac,ifelse(ctrl_label!=""," vs. ",""),ctrl_label))){
									beta_temp = beta_temp[which(pvalue[,5]<as.numeric(tclvalue(permutecut.val))),]
									peak_index_temp = peak_index[which(pvalue[,5]<as.numeric(tclvalue(permutecut.val))),]
									pvalue = read.table(sprintf("%s/Empirical p-value_%s%s%s.csv",outpath,fac,ifelse(ctrl_label!=""," vs. ",""),ctrl_label), sep = ",", header = T, quote = "\"")
									pfdr = pvalue
									
								

									png(filename = sprintf("%s/Volcanoplot_%s%s%s_epFDR_%d.png",outpath,fac,ifelse(ctrl_label!=""," vs. ",""),ctrl_label,top.num), width = 1440, height = 800)
									i=5
									par(mar = c(5,4,1,3))
									pch = rep(21, nrow(pfdr))
									pch[which(-log10(pfdr[,i])>=-log10(0.05))] = 24
									bg = rep("white", nrow(pfdr))
									bg[which(-log10(pfdr[,i])>=-log10(0.05))] = "green3"
									bg[order(-log10(pfdr[,i]), decreasing = T)[1:min(top.num, nrow(pfdr))]] = "red"
									plot(beta_temp[,4], -log10(pfdr[,i]), xlab = "", ylab = "", cex = 1, pch = pch, bg = bg)
									text(beta_temp[order(-log10(pfdr[,i]), decreasing = T)[1:min(top.num, nrow(pfdr))],4], -log10(pfdr[order(-log10(pfdr[,i]), decreasing = T)[1:min(top.num, nrow(pfdr))],i]), label = peak_index_temp[order(-log10(pfdr[,i]), decreasing = T)[1:min(top.num, nrow(pfdr))],1], cex = 1, pos = 1, adj = 0.2)

									abline(h = -log10(0.05), col = "green3")
									par(xpd = NA)
									text(par("usr")[1], -log10(0.05), col = "green3", label = substitute(paste(plain(-log[10]), "(0.05) (", test, " significant peaks)"), list(test = sum(pfdr[,i]<0.05))), font = 4, cex = 1.5, adj = c(0,-0.1))
									par(xpd = F)
									mtext("Fold Change (coefficient)", line = 2, side = 1, col="black", padj = 1, cex = 1.5)
									mtext(expression(-log[10](epFDR)), line = 3.5, side = 2, col="black", padj = 1, cex = 1.5)
									legend(x = "topleft", legend = c("Empirical p-FDR<0.05", paste("top", min(top.num, nrow(pfdr)), sep = "")), pch = 24, bty = "o", box.col = "black", cex = 1, pt.bg = c("green3", "red"), bg = rgb(1, 1, 1, alpha = 0.7))
									dev.off()
									
									png(filename = sprintf("%s/Volcanoplot_%s%s%s_epv_%d.png",outpath,fac,ifelse(ctrl_label!=""," vs. ",""),ctrl_label,top.num), width = 1440, height = 800)
									i=4
									par(mar = c(5,4,1,3))
									pch = rep(21, nrow(pvalue))
									pch[which(-log10(pvalue[,i])>=-log10(0.05))] = 24
									bg = rep("white", nrow(pvalue))
									bg[which(-log10(pvalue[,i])>=-log10(0.05))] = "green3"
									bg[order(-log10(pvalue[,i]), decreasing = T)[1:min(top.num, nrow(pvalue))]] = "red"
									plot(beta_temp[,4], -log10(pvalue[,i]), cex = 1, pch = pch, bg = bg, xlab = "", ylab = "")
									text(beta_temp[order(-log10(pvalue[,i]), decreasing = T)[1:min(top.num, nrow(pvalue))],4], -log10(pvalue[order(-log10(pvalue[,i]), decreasing = T)[1:min(top.num, nrow(pvalue))],i]), label = peak_index_temp[order(-log10(pvalue[,i]), decreasing = T)[1:min(top.num, nrow(pvalue))],1], cex = 1, pos = 1, adj = 0.2)
									abline(h = -log10(0.05), col = "green3")
									par(xpd = NA)
									text(par("usr")[1], -log10(0.05), col = "green3", label = substitute(paste(plain(-log[10]), "(0.05) (", test, " significant peaks)"), list(test = sum(pvalue[,i]<0.05))), font = 4, cex = 1.5, adj = c(0,-0.1))
									par(xpd = F)
									mtext("Fold Change (coefficient)", line = 2, side = 1, col="black", padj = 1, cex = 1.5)
									mtext(expression(-log[10](epv)), line = 3.5, side = 2, col="black", padj = 1, cex = 1.5)
									legend(x = "topleft", legend = c("Empirical p-value<0.05", paste("top", min(top.num, nrow(pvalue)), sep = "")), pch = 24, bty = "o", box.col = "black", cex = 1, pt.bg = c("green3", "red"), bg = rgb(1, 1, 1, alpha = 0.7))
									dev.off()

								}
							}
							tkmessageBox(title = "Analysis of Covariance", message = "Graphs are saved.", icon = "info", type = "ok")
							cat("ANCOVA - Graphs are saved.\n", sep = "")
							tkconfigure(tt, cursor = "arrow")
						})
						close.but <- tkbutton(fr_toppeak, text = "  Close  ", command = function(...){
							tkdestroy(toppeak)
							tkmessageBox(title="Analysis of Covariance", message = "ANCOVA is done.", icon = "info", type = "ok")
							cat("ANCOVA - ANCOVA is done.\n", sep = "")
						})
						
						tkgrid(plot.but, tklabel(fr_toppeak, text = "                              "), close.but)
						tkgrid(tklabel(toppeak, text = "", height = 0, font = fontIntro_para))
					}), side = "right")
					
					tkpack(tklabel(fr_next, text = "", height = 0, font = fontIntro_para), side = "bottom")
					tkconfigure(tt,cursor="arrow")
					tkmessageBox(title="Analysis of Covariance", message = "Please input the VIF upper bound for claiming a collinearity between batch effect and covariate(s).", icon = "info", type = "ok")
					cat("ANCOVA - Please input the VIF upper bound for claiming a collinearity between batch effect and covariate(s).\n", sep = "")

				}else{
					tkmessageBox(title = "Error", message = paste("Sample ID not match.\nPlease input the adequate ", paste(c("factor", "batch effect")[which(!c(all(!is.na(match(gsub("_[0-9]", "", colnames(peakabun)), Factor[,1]))), all(!is.na(match(colnames(peakabun), batch_effect[,1])))))], collapse = " and "), " file.", sep = ""), icon = "error", type = "ok")
					cat("ANCOVA(Error) - ", paste("Sample ID not match. Please input the adequate ", paste(c("factor", "batch effect")[which(!c(all(!is.na(match(gsub("_[0-9]", "", colnames(peakabun)), Factor[,1]))), all(!is.na(match(colnames(peakabun), batch_effect[,1])))))], collapse = " and "), " file", sep = ""), ".\n")
					tkfocus(dlg)
				}
			}else if((tclvalue(textbatchinput)=="") & !(tclvalue(textcova)=="")){  
				Covariate = read.table(tclvalue(textcova), header = T, quote = "\"", fill = T, sep = ifelse(grepl(".txt", tclvalue(textcova)), "\t", ","), as.is = T, na.string = c("NA", as.character(tclvalue(text_missing))))  
				cova_type = colnames(Covariate)
				cova_type = toupper(substring(cova_type[-1], 1, 1))
				
				peak_index = peakabun[,c(1:3)]
				peakabun = peakabun[,-c(1:3)]
				colnames(peakabun) = gsub("^\\d{4,8}_(.*?)|(^[A-Z]\\d{4,8})_(.*?)", "\\1", colnames(peakabun))
				Covariate[,1] = gsub("^\\d{4,8}_(.*?)|(^[A-Z]\\d{4,8})_(.*?)", "\\1", Covariate[,1])
				Factor[,1] = gsub("^\\d{4,8}_(.*?)|(^[A-Z]\\d{4,8})_(.*?)", "\\1", Factor[,1])
				
				if(all(!is.na(match(gsub("_[0-9]", "", colnames(peakabun)), Factor[,1]))) & all(!is.na(match(gsub("_[0-9]", "", colnames(peakabun)), Covariate[,1]))))
				{
					tkdestroy(dlg)
					tkconfigure(tt,cursor="watch")
					
					DRep = factor(gsub("(.*?)_(.*?)", "\\2", colnames(peakabun)))
					DInd <- gsub("_.*", "", colnames(peakabun))
					if(ctrl_label!=""){
						Factor = Factor[match(gsub("_[0-9]", "", colnames(peakabun)), Factor[,1]),2]
						Factor_label = levels(factor(Factor))
						Factor_label = Factor_label[-which(Factor_label==ctrl_label)]
					} else {
						Factor <- Factor[match(gsub("_[0-9]", "", colnames(peakabun)), Factor[,1]),]
						Factor <- subset(Factor,select=-1) 
						Factor_label <- colnames(Factor)
					}
					Covariate = if(ncol(Covariate)==2){
						as.data.frame(matrix(Covariate[match(gsub("_[0-9]", "", colnames(peakabun)), Covariate[,1]),-1], ncol = ncol(Covariate)-1, dimnames = list(NULL, colnames(Covariate)[2])))
					}else{
						as.data.frame(Covariate[match(gsub("_[0-9]", "", colnames(peakabun)), Covariate[,1]),-1])
					}
					
					for(x in 1:ncol(Covariate)){
						if(cova_type[x]=="D"){
							Covariate[,x] = factor(Covariate[,x])
						}
					}
					colnames(Covariate) = substring(colnames(Covariate), 2, nchar(colnames(Covariate)))
					
					
					
					  
					  
					peakabun = t(peakabun)
					
					if(ctrl_label!="")
						ancova_cov = as.vector(sapply(Factor_label, function(x) paste(colnames(Covariate), " (", x, " vs. ", ctrl_label, ")", sep = "")))
					else
						ancova_cov <- as.vector(sapply(Factor_label, function(x) paste(colnames(Covariate), " (", x, ")", sep = "")))
					ancova_var_idx <<- 1:length(ancova_cov)
					ancova_cov_idx <<- NULL
					cova_select <- tktoplevel(); if(isIcon) tk2ico.set(cova_select,icon)
					tkwm.title(cova_select, "Analysis of Covariance - Covariate Selection")
					choose_frame <- tkframe(cova_select)
					var_frame <- tkframe(choose_frame)
					var_scr <- tkscrollbar(var_frame, repeatinterval = 5, command = function(...) tkyview(var_tl,...))
					var_tl <- tklistbox(var_frame, height = 15, selectmode = "extended", yscrollcommand = function(...) tkset(var_scr,...), background = "white")
					tkgrid(tklabel(var_frame, text = "Variable"))
					tkgrid(var_tl, var_scr)
					tkgrid.configure(var_scr, rowspan = 4, sticky = "nsw")
					for(x in ancova_var_idx){
						tkinsert(var_tl, "end", ancova_cov[x])
					}
					cova_frame <- tkframe(choose_frame)
					cova_scr <- tkscrollbar(cova_frame, repeatinterval = 5, command = function(...) tkyview(cova_tl,...))
					cova_tl <- tklistbox(cova_frame, height = 15, selectmode = "extended", yscrollcommand = function(...) tkset(cova_scr, ...), background = "white")
					tkgrid(tklabel(cova_frame, text = "Covariate"))
					tkgrid(cova_tl, cova_scr)
					tkgrid.configure(cova_scr, rowspan = 4, sticky = "nsw")
					in_out <- tkframe(choose_frame)
					tkpack(tkbutton(in_out, text = "  >>  ", command = function(...){
						if(length(as.integer(tkcurselection(var_tl)))!=0){
							varIndex <- as.integer(tkcurselection(var_tl))
							for(x in 1:length(varIndex)){
								tkdelete(var_tl, varIndex[x]-x+1)
							}
							ancova_cov_idx <<- sort(c(ancova_cov_idx, ancova_var_idx[varIndex+1]))
							for(x in varIndex){
								tkinsert(cova_tl, which(ancova_cov_idx==ancova_var_idx[x+1])-1, ancova_cov[ancova_var_idx[x+1]])
							}
							ancova_var_idx <<- ancova_var_idx[-(varIndex+1)]
						}
						tkconfigure(but.covSel,state=ifelse(length(ancova_cov_idx)>0,"normal","disable"))
					}))
					tkpack(tklabel(in_out, text = "     "))
					tkpack(tkbutton(in_out, text = "  <<  ", command = function(...){
						if(length(as.integer(tkcurselection(cova_tl)))!=0){
							covIndex <- as.integer(tkcurselection(cova_tl))
							for(x in 1:length(covIndex)){
								tkdelete(cova_tl, covIndex[x]-x+1)
							}
							ancova_var_idx <<- sort(c(ancova_var_idx, ancova_cov_idx[covIndex+1]))
							for(x in covIndex){
								tkinsert(var_tl, which(ancova_var_idx==ancova_cov_idx[x+1])-1, ancova_cov[ancova_cov_idx[x+1]])
							}
							ancova_cov_idx <<- ancova_cov_idx[-(covIndex+1)]
						}
						tkconfigure(but.covSel,state=ifelse(length(ancova_cov_idx)>0,"normal","disable"))
					}))

					tkgrid(var_frame, in_out, cova_frame, padx = 10)
					tkgrid(tklabel(choose_frame, text = "", height = 0, font = fontIntro_para))
					tkgrid(choose_frame)
					
					
					fr_next <- tkframe(cova_select)
					tkpack(but.covSel <- tkbutton(fr_next, text = "  Next  ", state="disabled", command = function(...){

						if(check) tkmessageBox(title="Analysis of Covariance", message = sprintf("Instead of nested ANCOVA, ANCOVA is applied\nbecause %d samples have only one replicate.",check), icon = "info", type = "ok")
				
						tkdestroy(cova_select)
						tkconfigure(tt,cursor="watch")
						
						if(sum(ancova_cov_idx)){
							ancova_cov[-ancova_cov_idx] = 0
							ancova_cov[ancova_cov_idx] = 1
						}
						ancova_cov = matrix(as.numeric(ancova_cov), nrow = length(colnames(Covariate)), dimnames = list(colnames(Covariate), Factor_label))
						nested_ancova <- tktoplevel(); if(isIcon) tk2ico.set(nested_ancova,icon)
						tkwm.title(nested_ancova, "Analysis of Covariance - Model")
						tkpack(fr_head <- tkframe(nested_ancova), side = "top", fill = "x")
						tkpack(tklabel(fr_head, text = "Please input an ANCOVA model."), side = "left")
						tkpack(fr_input <- tkframe(nested_ancova), side = "top", fill = "x")
						for(grp in 1:ncol(ancova_cov)){
							tkpack(fr_input.1 <- tkframe(nested_ancova), side = "top", fill = "x")
							tkpack(fr_input.2 <- tkframe(nested_ancova), side = "top", fill = "x")
							if(all(ancova_cov[,grp]==0)){
								tkpack(tklabel(fr_input.1, text = paste("     Candidate covariates:  ", paste(ifelse(ctrl_label!="","Group , ","Quan , "),ifelse(check,"","Rep , "), paste(rownames(ancova_cov)[which(ancova_cov[,Factor_label[grp]]==1)], collapse = " , "), sep = ""), sep = "")), side = "left")
								assign(paste("model_", Factor_label[grp], sep = ""), tclVar(paste(c(if(check) NULL else ifelse(ctrl_label!="","Group:Rep","Quan:Rep"), rownames(ancova_cov)[which(ancova_cov[,Factor_label[grp]]==1)]), collapse = " + ")))
								if(ctrl_label!="") tkpack(tklabel(fr_input.2, text = paste("     ", Factor_label[grp], " vs. ", ctrl_label, ":   Group + ", sep = "")), tkentry(fr_input.2, width = "80", textvariable = tclVar("Group:Rep"), state = "disable"), tklabel(fr_input.2, text = "     "), side = "left")
								else tkpack(tklabel(fr_input.2, text = paste("     ", Factor_label[grp], " vs. ", ctrl_label, ":   Quan + ", sep = "")), tkentry(fr_input.2, width = "80", textvariable = tclVar("Quan:Rep"), state = "disable"), tklabel(fr_input.2, text = "     "), side = "left")

							}else{
								tkpack(tklabel(fr_input.1, text = paste("Candidate covariates:  ", paste(ifelse(ctrl_label!="","Group , ","Quan , "),ifelse(check,"","Rep , "), paste(rownames(ancova_cov)[which(ancova_cov[,Factor_label[grp]]==1)], collapse = " , "), sep = ""), sep = ""), wraplength = 750, justify = "left", padx = 17), side = "left")
								assign(paste("model_", Factor_label[grp], sep = ""), tclVar(paste(c(if(check) NULL else ifelse(ctrl_label!="","Group:Rep","Quan:Rep"), rownames(ancova_cov)[which(ancova_cov[,Factor_label[grp]]==1)]), collapse = " + ")))
								if(ctrl_label!="") tkpack(tklabel(fr_input.2, text = paste("     ", Factor_label[grp], " vs. ", ctrl_label, ":   Group + ", sep = "")), tkentry(fr_input.2, width = "80", textvariable = get(paste("model_", Factor_label[grp], sep = "")), bg = "white"), tklabel(fr_input.2, text = "     "), side = "left")
								else tkpack(tklabel(fr_input.2, text = paste("     ", Factor_label[grp], ":   Quan + ", sep = "")), tkentry(fr_input.2, width = "80", textvariable = get(paste("model_", Factor_label[grp], sep = "")), bg = "white"), tklabel(fr_input.2, text = "     "), side = "left")
							}
						}
						tkpack(tklabel(nested_ancova, text = "", height = 0, font = fontIntro_para), side = "top")
						tkpack(fr_note <- tkframe(nested_ancova), side = "top", fill = "x")
						tkpack(tklabel(fr_note, text = "", height = 0, font = fontIntro_para), side = "bottom")
						tkpack(tklabel(fr_note, text = "     Note: Notation ''A:B'' indicates that variable B nested in variable A.", justify = "left"), side = "left")
						tkpack(tklabel(fr_note, text = "     "), side = "right")
						tkpack(tkbutton(fr_note, text = "  Next  ", command = function(...){
							tkconfigure(tt,cursor="watch")
							model = rep(NA, ncol(ancova_cov))
							for(grp in 1:ncol(ancova_cov)){
								if(all(ancova_cov[,grp]==0)) next
								model[grp] = gsub("\n", "", tclvalue(get(paste("model_", Factor_label[grp], sep = ""))))
							}
							tkdestroy(nested_ancova)

							op <- options(contrasts = c("contr.helmert", "contr.poly"))
							pb <- tkProgressBar("ANCOVA - Please wait for ANCOVA processing", "0% done", 0, 100, 0, width = 500)
							cat("ANCOVA - Please wait for ANCOVA processing 0 ")
							if(tclvalue(permute.val)=="1" & ncpu>1) cl = makeCluster(ncpu, type = "SOCK")
							for(fa in 1:length(Factor_label)){
								fac = Factor_label[fa]
								check_model = any(ancova_cov[,fac]==1)
								model_choose = model[fa]
								if(ctrl_label!=""){
									group_idx = which(Factor %in% c(fac, ctrl_label))
									Factor_temp = factor(Factor[group_idx])
								} else {
									group_idx <- 1:nrow(Factor)
									Factor_temp <- Factor[,fa]
								}
								Rep = factor(DRep[group_idx])
								Ind = factor(DInd[group_idx],level=unique(DInd[group_idx]))

								
								cov_idx = which(ancova_cov[,fac]==1)
								
								if(check_model){
									temp_all = as.data.frame(Covariate[group_idx, cov_idx])
									colnames(temp_all) = colnames(Covariate)[cov_idx]

									options(contrasts = c("contr.helmert", "contr.poly"))
									anova <- apply(peakabun[group_idx,],2,function(x) {
										nLevel <- apply(cbind(Factor_temp,temp_all)[!is.na(x),],2,function(y) length(unique(y)))
										if(any(nLevel<2)) NA else aov(as.formula(paste("x ~ Factor_temp + ", gsub("Group|Quan", "Factor_temp", model[fa]))), data = temp_all)
									})
									pvalue <- do.call("rbind",sapply(anova,function(x) if(all(is.na(x))) NA else if(summary(x)[[1]]["Residuals","Sum Sq"]<2.2e-16) NA else Anova(x,type="III",singular.ok=T)[,"Pr(>F)"],simplify=F))
									colnames(pvalue) <- gsub(" ","",gsub("_temp","",gsub("Factor_temp", sprintf("%s%s%s",fac,ifelse(ctrl_label!="","_",""),ctrl_label), c("(Intercept)",rownames(summary(anova[[which(!apply(pvalue,1,function(x) all(is.na(x))))[1]]])[[1]])))))
									options(contrasts = c("contr.SAS", "contr.poly"))
									if(ctrl_label!="") contrasts(Factor_temp) = contrasts(Factor_temp)[c(which(levels(Factor_temp)!=ctrl_label), which(levels(Factor_temp)==ctrl_label)),]
									beta_temp <- apply(peakabun[group_idx,],2,function(x) {
										nLevel <- apply(cbind(Factor_temp,temp_all)[!is.na(x),],2,function(y) length(unique(y)))
										if(any(nLevel<2)) NA else lm(as.formula(paste("x ~ Factor_temp + ", gsub("Group|Quan", "Factor_temp", model[fa]))), data = temp_all)$coefficients
									})
									beta_temp = switch(class(beta_temp), "matrix" = t(beta_temp), "list"=do.call(rbind, beta_temp))
									resi_temp = apply(peakabun[group_idx,], 2, function(x){
										nLevel <- apply(cbind(Factor_temp,temp_all)[!is.na(x),],2,function(y) length(unique(y)))
										temp <- rep(NA, length(x))
										if(all(nLevel>=2)) temp[!is.na(x)&rowSums(is.na(temp_all))==0] <- lm(as.formula(paste("x ~ ", gsub("Group|Quan", "Factor_temp", model[fa]))), data = temp_all)$residuals
										return(temp)
									})
									resi_temp = switch(class(resi_temp), "matrix" = t(resi_temp), "list"=do.call(rbind, resi_temp))
							
								}else{

									options(contrasts = c("contr.helmert", "contr.poly"))
									anova <- apply(peakabun[group_idx,],2,function(x) {
										nLevel <- apply(cbind(Factor_temp,temp_all)[!is.na(x),],2,function(y) length(unique(y)))
										if(any(nLevel<2)) NA else aov(as.formula(ifelse(check,"x ~ Factor_temp","x ~ Factor_temp + Factor_temp:Rep")))
									})
									pvalue <- do.call("rbind",sapply(anova,function(x) if(all(is.na(x))) NA else if(summary(x)[[1]]["Residuals","Sum Sq"]<2.2e-16) NA else Anova(x,type="III",singular.ok=T)[,"Pr(>F)"],simplify=F))
									colnames(pvalue) <- gsub(" ","",gsub("_temp","",gsub("Factor_temp", sprintf("%s%s%s",fac,ifelse(ctrl_label!="","_",""),ctrl_label), c("(Intercept)",rownames(summary(anova[[which(!apply(pvalue,1,function(x) all(is.na(x))))[1]]])[[1]])))))
									options(contrasts = c("contr.SAS", "contr.poly"))
									if(ctrl_label!="") contrasts(Factor_temp) = contrasts(Factor_temp)[c(which(levels(Factor_temp)!=ctrl_label), which(levels(Factor_temp)==ctrl_label)),]
									beta_temp <- apply(peakabun[group_idx,],2,function(x) {
										nLevel <- length(unique(Factor_temp[!is.na(x)]))
										if(any(nLevel<2)) NA else lm(as.formula(ifelse(check,"x ~ Factor_temp","x ~ Factor_temp + Factor_temp:Rep")))$coefficients
									})
									beta_temp = switch(class(beta_temp), "matrix" = t(beta_temp), "list"=do.call(rbind, beta_temp))
									resi_temp = apply(peakabun[group_idx,], 2, function(x){
										nLevel <- apply(cbind(Factor_temp,temp_all)[!is.na(x),],2,function(y) length(unique(y)))
										temp <- rep(NA, length(x))
										if(all(nLevel>=2)) temp[!is.na(x)&rowSums(is.na(temp_all))==0] <- lm(as.formula(ifelse(check,"x ~ 1","x ~ Factor_temp:Rep")))$residuals
										return(temp)
									})
									resi_temp = switch(class(resi_temp), "matrix" = t(resi_temp), "list"=do.call(rbind, resi_temp))

								}
								pvalue = pvalue[,-ncol(pvalue)]
								write.table(cbind(peak_index, pvalue), file = sprintf("%s/p-value_%s%s%s_details.csv",outpath,fac,ifelse(ctrl_label!=""," vs. ",""),ctrl_label), col.names = T, row.names = F, quote = 1, sep = ",")
								pvalue = as.data.frame(pvalue[,2])
								beta_temp = as.data.frame(beta_temp[,2])
								resi_temp = as.data.frame(resi_temp); colnames(resi_temp) <- rownames(peakabun[group_idx,])
								colnames(pvalue) = sprintf("%s%s%s",fac,ifelse(ctrl_label!="","_",""),ctrl_label)
								pfdr = apply(pvalue, 2, function(x) p.adjust(x, method = "BH"))
								colnames(pvalue) = paste(colnames(pvalue), "(pval)", sep = "")
								colnames(pfdr) = paste(colnames(pfdr), "(pfdr)", sep = "")
								colnames(beta_temp) = sprintf("%s%s%s(coef)",fac,ifelse(ctrl_label!="","_",""),ctrl_label)
								
								write.table(cbind(peak_index, pfdr), file = sprintf("%s/pFDR_%s%s%s.csv",outpath,fac,ifelse(ctrl_label!=""," vs. ",""),ctrl_label), col.names = T, row.names = F, quote = 1, sep = ",")
								write.table(cbind(peak_index, beta_temp), file = sprintf("%s/coef_%s%s%s.csv",outpath,fac,ifelse(ctrl_label!=""," vs. ",""),ctrl_label), col.names = T, row.names = F, quote = 1, sep = ",")
								write.table(cbind(peak_index, resi_temp), file = sprintf("%s/resi_%s%s%s_forPLS.csv",outpath,fac,ifelse(ctrl_label!=""," vs. ",""),ctrl_label), col.names = T, row.names = F, quote = 1, sep = ",")
								if(tclvalue(permute.val)!="1") write.table(cbind(peak_index, resi_temp)[pfdr<.05,], file = sprintf("%s/resi_pFDRsig_%s%s%s_forPLS.csv",outpath,fac,ifelse(ctrl_label!=""," vs. ",""),ctrl_label), col.names = T, row.names = F, quote = 1, sep = ",")
								
								if(tclvalue(permute.val)=="1"){
									if(sum(pvalue[,1]<as.numeric(tclvalue(permutecut.val)),na.rm=T)>0){
										peak_sig = which(pvalue[,1]<as.numeric(tclvalue(permutecut.val)))
										pvalue_sig = pvalue[peak_sig, 1]
										peakabun_sig = as.data.frame(peakabun[group_idx,peak_sig])
										if(ncpu==1){
											pvalue_permute = sapply(1:as.numeric(tclvalue(permute.times)), function(perm){
												Factor_temp_perm = sample(Factor_temp)
												
												options(contrasts = c("contr.helmert", "contr.poly"))
												if(check_model){
													temp_all2 <- temp_all
													
														
													
														
														

														
														
														
														
													
													anova_temp <- apply(peakabun_sig,2,function(x) {
														nLevel <- length(unique(Factor_temp_perm[!is.na(x)]))
														if(any(nLevel<2)) NA else aov(as.formula(paste("x ~ Factor_temp_perm + ", gsub("Group", "Factor_temp_perm", model_choose))), data = temp_all2)
													})
													pvalue_temp <- do.call("cbind",sapply(anova_temp,function(x) if(all(is.na(x))) NA else if(summary(x)[[1]]["Residuals","Sum Sq"]<2.2e-16) NA else Anova(x,type="III",singular.ok=T)[,"Pr(>F)"],simplify=F))
												}else{
													
														
													
														
														
														
													
													anova_temp <- apply(peakabun_sig,2,function(x) {
														nLevel <- length(unique(Factor_temp_perm[!is.na(x)]))
														if(any(nLevel<2)) NA else aov(as.formula(ifelse(check,"x ~ Factor_temp_perm","x ~ Factor_temp_perm + Factor_temp_perm:Rep")))
														
													})
													pvalue_temp <- do.call("cbind",sapply(anova_temp,function(x) if(all(is.na(x))) NA else if(summary(x)[[1]]["Residuals","Sum Sq"]<2.2e-16) NA else Anova(x,type="III",singular.ok=T)[,"Pr(>F)"],simplify=F))
												}
												info <- sprintf("%d%%", round(100*fa/length(Factor_label)*(perm/as.numeric(tclvalue(permute.times)))))
												setTkProgressBar(pb, value = round(100*fa/length(Factor_label)*(perm/as.numeric(tclvalue(permute.times)))), sprintf("ANCOVA-Please wait for ANCOVA processing (%s done)", info), paste(fac, "-", ctrl_label, info, sep = " "))
												return(pvalue_temp[2,]<pvalue_sig)
											})
										}else if(ncpu>1){
											if(check_model){
												clusterExport(cl, c("Factor_temp", "Ind", "Rep", "peakabun_sig", "pvalue_sig", "temp_all", "check_model", "Anova", "model_choose", "check", "osIsWin"), envir=environment())
												
											}else{
												clusterExport(cl, c("Factor_temp", "Ind", "Rep", "peakabun_sig", "pvalue_sig", "check_model", "Anova", "model_choose", "check", "osIsWin"), envir=environment())
												
											}
											pvalue_permute = parSapply(cl, 1:as.numeric(tclvalue(permute.times)), function(perm){

												pb <- tkProgressBar(sprintf("ANCOVA-Parallel task (permutation %d)",perm), "0% done", 0, 100, 0, width = 500)
													Factor_temp_perm = sample(Factor_temp)
													
													options(contrasts = c("contr.helmert", "contr.poly"))
													if(check_model){
														temp_all2 <- temp_all
														
															
														
															
															

															
															
															
															
														
														anova_temp <- apply(peakabun_sig,2,function(x) {
															nLevel <- length(unique(Factor_temp_perm[!is.na(x)]))
															if(any(nLevel<2)) NA else aov(as.formula(paste("x ~ Factor_temp_perm + ", gsub("Group", "Factor_temp_perm", model_choose))), data = temp_all2)
														})
														pvalue_temp <- do.call("cbind",sapply(anova_temp,function(x) if(all(is.na(x))) NA else if(summary(x)[[1]]["Residuals","Sum Sq"]<2.2e-16) NA else Anova(x,type="III",singular.ok=T)[,"Pr(>F)"],simplify=F))
													}else{
														
															
														
															
															
															
														
														anova_temp <- apply(peakabun_sig,2,function(x) {
															nLevel <- length(unique(Factor_temp_perm[!is.na(x)]))
															if(any(nLevel<2)) NA else aov(as.formula(ifelse(check,"x ~ Factor_temp_perm","x ~ Factor_temp_perm + Factor_temp_perm:Rep")))
															
														})
														pvalue_temp <- do.call("cbind",sapply(anova_temp,function(x) if(all(is.na(x))) NA else if(summary(x)[[1]]["Residuals","Sum Sq"]<2.2e-16) NA else Anova(x,type="III",singular.ok=T)[,"Pr(>F)"],simplify=F))
													}
												if(osIsWin) for(i in 1:100) setTkProgressBar(pb, value = i, sprintf("ANCOVA-Parallel task (permutation %d)", perm), "")
												return(pvalue_temp[2,]<pvalue_sig)
											})
										}
										empvalue = as.data.frame((rowSums(pvalue_permute)+1)/(as.numeric(tclvalue(permute.times))+1))
										empvalue = cbind(empvalue, p.adjust(c(empvalue[,1], rep(max(c(empvalue[,1], as.numeric(tclvalue(permutecut.val)))), nrow(peak_index)-nrow(empvalue))), method = "BH")[1:nrow(empvalue)])
										colnames(empvalue) = c(sprintf("%s%s%s(epv)",fac,ifelse(ctrl_label!="","_",""),ctrl_label), sprintf("%s%s%s(epFDR)",fac,ifelse(ctrl_label!="","_",""),ctrl_label))
										write.table(cbind(peak_index[peak_sig,], empvalue), file = sprintf("%s/Empirical p-value_%s%s%s.csv",outpath,fac,ifelse(ctrl_label!=""," vs. ",""),ctrl_label), col.names = T, row.names = F, quote = 1, sep = ",")
										write.table(cbind(peak_index, resi_temp)[peak_sig,][empvalue[,2]<.05,], file = sprintf("%s/resi_EmpiricalpFDRsig_%s%s%s_forPLS.csv",outpath,fac,ifelse(ctrl_label!=""," vs. ",""),ctrl_label), col.names = T, row.names = F, quote = 1, sep = ",")
									}else{
										tkmessageBox(title = "Warning", message = "No significant peaks", icon = "warning", type = "ok")
									}
								}
								
								info <- sprintf("%d%%", round(100*fa/length(Factor_label)))
								setTkProgressBar(pb, value = round(100*fa/length(Factor_label)), sprintf("ANCOVA - Please wait for ANCOVA processing (%s done)", info), paste(fac, "-", ctrl_label, info, sep = " "))
								if(round(100*fa/length(Factor_label))<100){
									cat(round(100*fa/length(Factor_label)), " ", sep = "")
								}else{
									cat(round(100*fa/length(Factor_label)), " \n", sep = "")
								}
								Sys.sleep(0.1)
							}
							if(tclvalue(permute.val)=="1" & ncpu>1) stopCluster(cl)
							options(op)
							ancova_cov = cbind(c("Cov", rownames(ancova_cov)), rbind(colnames(ancova_cov), ancova_cov))
							write.table(ancova_cov, file = paste(outpath, "/covariate_select.csv", sep = ""), sep = ",", col.names = F, row.names = F, quote = F)
							setTkProgressBar(pb, value = 100, "ANCOVA - Please wait for ANCOVA processing (100% done)", "Finished 100%")
							Sys.sleep(1)
							close(pb)
							tkconfigure(tt,cursor="arrow")
							
							toppeak = tktoplevel(); if(isIcon) tk2ico.set(toppeak,icon)
							tkwm.title(toppeak, "Top peaks")
							fr_top = tkframe(toppeak)
							tkgrid(fr_top)
							top.lab <- tklabel(fr_top, text = "    How many top peaks do you want to paint red?  ")
							top.num = tclVar("50")
							top.entry <- tkentry(fr_top, width = 6, textvariable = top.num, bg = "white",validatecommand="string is integer %P", validate = "all")
							
							tkgrid(top.lab, top.entry, tklabel(fr_top, text = paste("(1 ~ ", ncol(peakabun), ")     ", sep = "")))
							tkgrid(tklabel(toppeak, text = "", height = 0, font = fontIntro_para))
							fr_toppeak = tkframe(toppeak)
							tkgrid(fr_toppeak)
							plot.but <- tkbutton(fr_toppeak, text = "  Plot  ", command = function(...){
								tkconfigure(tt,cursor="watch")
								top.num = as.numeric(tclvalue(top.num))
								if(!top.num%in%seq.int(1, ncol(peakabun))) top.num = 50
								for(fa in 1:length(Factor_label)){
									fac = Factor_label[fa]
									pfdr = read.table(sprintf("%s/pFDR_%s%s%s.csv",outpath,fac,ifelse(ctrl_label!=""," vs. ",""),ctrl_label), sep = ",", header = T, quote = "\"")
									pvalue = read.table(sprintf("%s/p-value_%s%s%s_details.csv",outpath,fac,ifelse(ctrl_label!=""," vs. ",""),ctrl_label), sep = ",", header = T, quote = "\"")
									beta_temp = read.table(sprintf("%s/coef_%s%s%s.csv",outpath,fac,ifelse(ctrl_label!=""," vs. ",""),ctrl_label), sep = ",", header = T, quote = "\"")
									
								

									png(filename = sprintf("%s/Volcanoplot_%s%s%s_pFDR_%d.png",outpath,fac,ifelse(ctrl_label!=""," vs. ",""),ctrl_label,top.num), width = 1440, height = 800)
									i=4
									par(mar = c(5,4,1,3))
									pch = rep(21, nrow(pfdr))
									pch[which(-log10(pfdr[,i])>=-log10(0.05))] = 24
									bg = rep("white", nrow(pfdr))
									bg[which(-log10(pfdr[,i])>=-log10(0.05))] = "green3"
									bg[order(-log10(pfdr[,i]), decreasing = T)[1:top.num]] = "red"
									plot(beta_temp[,4], -log10(pfdr[,i]), xlab = "", ylab = "", cex = 1, pch = pch, bg = bg)
									text(beta_temp[order(-log10(pfdr[,i]), decreasing = T)[1:top.num],4], -log10(pfdr[order(-log10(pfdr[,i]), decreasing = T)[1:top.num],i]), label = peak_index[order(-log10(pfdr[,i]), decreasing = T)[1:top.num],1], cex = 1, pos = 1, adj = 0.2)

									abline(h = -log10(0.05), col = "green3")
									par(xpd = NA)
									text(par("usr")[1], -log10(0.05), col = "green3", label = substitute(paste(plain(-log[10]), "(0.05) (", test, " significant peaks)"), list(test = sum(pfdr[,i]<0.05))), font = 4, cex = 1.5, adj = c(0,-0.1))
									par(xpd = F)
									mtext("Fold Change (coefficient)", line = 2, side = 1, col="black", padj = 1, cex = 1.5)
									mtext(expression(-log[10](pFDR)), line = 3.5, side = 2, col="black", padj = 1, cex = 1.5)
									legend(x = "topleft", legend = c("p-FDR<0.05", paste("top", top.num, sep = "")), pch = 24, bty = "o", box.col = "black", cex = 1, pt.bg = c("green3", "red"), bg = rgb(1, 1, 1, alpha = 0.7))
									dev.off()
									
									png(filename = sprintf("%s/Volcanoplot_%s%s%s_pv_%d.png",outpath,fac,ifelse(ctrl_label!=""," vs. ",""),ctrl_label,top.num), width = 1440, height = 800)
									i=5
									par(mar = c(5,4,1,3))
									pch = rep(21, nrow(pvalue))
									pch[which(-log10(pvalue[,i])>=-log10(0.05))] = 24
									bg = rep("white", nrow(pvalue))
									bg[which(-log10(pvalue[,i])>=-log10(0.05))] = "green3"
									bg[order(-log10(pvalue[,i]), decreasing = T)[1:top.num]] = "red"
									plot(beta_temp[,4], -log10(pvalue[,i]), cex = 1, pch = pch, bg = bg, xlab = "", ylab = "")
									text(beta_temp[order(-log10(pvalue[,i]), decreasing = T)[1:top.num],4], -log10(pvalue[order(-log10(pvalue[,i]), decreasing = T)[1:top.num],i]), label = peak_index[order(-log10(pvalue[,i]), decreasing = T)[1:top.num],1], cex = 1, pos = 1, adj = 0.2)
									abline(h = -log10(0.05), col = "green3")
									par(xpd = NA)
									text(par("usr")[1], -log10(0.05), col = "green3", label = substitute(paste(plain(-log[10]), "(0.05) (", test, " significant peaks)"), list(test = sum(pvalue[,i]<0.05))), font = 4, cex = 1.5, adj = c(0,-0.1))
									par(xpd = F)
									mtext("Fold Change (coefficient)", line = 2, side = 1, col="black", padj = 1, cex = 1.5)
									mtext(expression(-log[10](pv)), line = 3.5, side = 2, col="black", padj = 1, cex = 1.5)
									legend(x = "topleft", legend = c("p-value<0.05", paste("top", top.num, sep = "")), pch = 24, bty = "o", box.col = "black", cex = 1, pt.bg = c("green3", "red"), bg = rgb(1, 1, 1, alpha = 0.7))
									dev.off()

									
									if(tclvalue(permute.val)=="1" & file.exists(sprintf("%s/Empirical p-value_%s%s%s.csv",outpath,fac,ifelse(ctrl_label!=""," vs. ",""),ctrl_label))){
										beta_temp = beta_temp[which(pvalue[,5]<as.numeric(tclvalue(permutecut.val))),]
										peak_index_temp = peak_index[which(pvalue[,5]<as.numeric(tclvalue(permutecut.val))),]
										pvalue = read.table(sprintf("%s/Empirical p-value_%s%s%s.csv",outpath,fac,ifelse(ctrl_label!=""," vs. ",""),ctrl_label), sep = ",", header = T, quote = "\"")
										pfdr = pvalue
										
									

										png(filename = sprintf("%s/Volcanoplot_%s%s%s_epFDR_%d.png",outpath,fac,ifelse(ctrl_label!=""," vs. ",""),ctrl_label,top.num), width = 1440, height = 800)
										i=5
										par(mar = c(5,4,1,3))
										pch = rep(21, nrow(pfdr))
										pch[which(-log10(pfdr[,i])>=-log10(0.05))] = 24
										bg = rep("white", nrow(pfdr))
										bg[which(-log10(pfdr[,i])>=-log10(0.05))] = "green3"
										bg[order(-log10(pfdr[,i]), decreasing = T)[1:min(top.num, nrow(pfdr))]] = "red"
										plot(beta_temp[,4], -log10(pfdr[,i]), xlab = "", ylab = "", cex = 1, pch = pch, bg = bg)
										text(beta_temp[order(-log10(pfdr[,i]), decreasing = T)[1:min(top.num, nrow(pfdr))],4], -log10(pfdr[order(-log10(pfdr[,i]), decreasing = T)[1:min(top.num, nrow(pfdr))],i]), label = peak_index_temp[order(-log10(pfdr[,i]), decreasing = T)[1:min(top.num, nrow(pfdr))],1], cex = 1, pos = 1, adj = 0.2)

										abline(h = -log10(0.05), col = "green3")
										par(xpd = NA)
										text(par("usr")[1], -log10(0.05), col = "green3", label = substitute(paste(plain(-log[10]), "(0.05) (", test, " significant peaks)"), list(test = sum(pfdr[,i]<0.05))), font = 4, cex = 1.5, adj = c(0,-0.1))
										par(xpd = F)
										mtext("Fold Change (coefficient)", line = 2, side = 1, col="black", padj = 1, cex = 1.5)
										mtext(expression(-log[10](epFDR)), line = 3.5, side = 2, col="black", padj = 1, cex = 1.5)
										legend(x = "topleft", legend = c("Empirical p-FDR<0.05", paste("top", min(top.num, nrow(pfdr)), sep = "")), pch = 24, bty = "o", box.col = "black", cex = 1, pt.bg = c("green3", "red"), bg = rgb(1, 1, 1, alpha = 0.7))
										dev.off()
										
										png(filename = sprintf("%s/Volcanoplot_%s%s%s_epv_%d.png",outpath,fac,ifelse(ctrl_label!=""," vs. ",""),ctrl_label,top.num), width = 1440, height = 800)
										i=4
										par(mar = c(5,4,1,3))
										pch = rep(21, nrow(pvalue))
										pch[which(-log10(pvalue[,i])>=-log10(0.05))] = 24
										bg = rep("white", nrow(pvalue))
										bg[which(-log10(pvalue[,i])>=-log10(0.05))] = "green3"
										bg[order(-log10(pvalue[,i]), decreasing = T)[1:min(top.num, nrow(pvalue))]] = "red"
										plot(beta_temp[,4], -log10(pvalue[,i]), cex = 1, pch = pch, bg = bg, xlab = "", ylab = "")
										text(beta_temp[order(-log10(pvalue[,i]), decreasing = T)[1:min(top.num, nrow(pvalue))],4], -log10(pvalue[order(-log10(pvalue[,i]), decreasing = T)[1:min(top.num, nrow(pvalue))],i]), label = peak_index_temp[order(-log10(pvalue[,i]), decreasing = T)[1:min(top.num, nrow(pvalue))],1], cex = 1, pos = 1, adj = 0.2)
										abline(h = -log10(0.05), col = "green3")
										par(xpd = NA)
										text(par("usr")[1], -log10(0.05), col = "green3", label = substitute(paste(plain(-log[10]), "(0.05) (", test, " significant peaks)"), list(test = sum(pvalue[,i]<0.05))), font = 4, cex = 1.5, adj = c(0,-0.1))
										par(xpd = F)
										mtext("Fold Change (coefficient)", line = 2, side = 1, col="black", padj = 1, cex = 1.5)
										mtext(expression(-log[10](epv)), line = 3.5, side = 2, col="black", padj = 1, cex = 1.5)
										legend(x = "topleft", legend = c("Empirical p-value<0.05", paste("top", min(top.num, nrow(pvalue)), sep = "")), pch = 24, bty = "o", box.col = "black", cex = 1, pt.bg = c("green3", "red"), bg = rgb(1, 1, 1, alpha = 0.7))
										dev.off()

									}
								}
								tkmessageBox(title = "Analysis of Covariance", message = "Graphs are saved.", icon = "info", type = "ok")
								cat("ANCOVA - Graphs are saved.\n", sep = "")
								tkconfigure(tt,cursor="arrow")
							})
							close.but <- tkbutton(fr_toppeak, text = "  Close  ", command = function(...){
								tkdestroy(toppeak)
								tkmessageBox(title="Analysis of Covariance", message = "ANCOVA is done.", icon = "info", type = "ok")
								cat("ANCOVA - ANCOVA is done.\n", sep = "")
							})
							
							
							tkgrid(plot.but, tklabel(fr_toppeak, text = "                              "), close.but)
							tkgrid(tklabel(toppeak, text = "", height = 0, font = fontIntro_para))
						}), side = "right", fill = "y")
						tkconfigure(tt,cursor="arrow")
					}))
					tkgrid(fr_next, sticky = "s", padx = 10)
					tkgrid(tklabel(cova_select, text = "", height = 0, font = fontIntro_para))
					tkmessageBox(title = "ANCOVA - Covariates Select", message = "Please select covariates in ANCOVA.", icon = "info", type = "ok")
					cat("ANCOVA - Please select covariates in ANCOVA.\n")
				}else{
					tkmessageBox(title = "Error", message = paste("Sample ID not match.\nPlease input the adequate ", paste(c("factor", "covariates")[which(!c(all(!is.na(match(gsub("_[0-9]", "", colnames(peakabun)), Factor[,1]))), all(!is.na(match(gsub("_[0-9]", "", colnames(peakabun)), Covariate[,1])))))], collapse = " and "), " file.", sep = ""), icon = "error", type = "ok")
					cat("ANCOVA(Error) - ", paste("Sample ID not match. Please input the adequate ", paste(c("factor", "covariates")[which(!c(all(!is.na(match(gsub("_[0-9]", "", colnames(peakabun)), Factor[,1]))), all(!is.na(match(gsub("_[0-9]", "", colnames(peakabun)), Covariate[,1])))))], collapse = " and "), " file", sep = ""), ".\n")
					tkfocus(dlg)
				}
			}else{  
				Covariate = read.table(tclvalue(textcova), header = T, quote = "\"", fill = T, sep = ifelse(grepl(".txt", tclvalue(textcova)), "\t", ","), as.is = T, na.string = c("NA", as.character(tclvalue(text_missing))))  
				cova_type = colnames(Covariate)
				cova_type = toupper(substring(cova_type[-1], 1, 1))
				batch_effect = read.table(tclvalue(textbatchinput), header = T, quote = "\"", sep = ifelse(grepl(".txt", tclvalue(textbatchinput)), "\t", ","), as.is = T, na.string = c("NA", as.character(tclvalue(text_missing))))

				peak_index = peakabun[,c(1:3)]
				peakabun = peakabun[,-c(1:3)]
				colnames(peakabun) = gsub("^\\d{4,8}_(.*?)|(^[A-Z]\\d{4,8})_(.*?)", "\\1", colnames(peakabun))
				Covariate[,1] = gsub("^\\d{4,8}_(.*?)|(^[A-Z]\\d{4,8})_(.*?)", "\\1", Covariate[,1])
				Factor[,1] = gsub("^\\d{4,8}_(.*?)|(^[A-Z]\\d{4,8})_(.*?)", "\\1", Factor[,1])
				batch_effect[,1] = gsub("^\\d{4,8}_(.*?)|(^[A-Z]\\d{4,8})_(.*?)", "\\1", batch_effect[,1])
				
				if(all(!is.na(match(gsub("_[0-9]", "", colnames(peakabun)), Factor[,1]))) & all(!is.na(match(gsub("_[0-9]", "", colnames(peakabun)), Covariate[,1]))) & all(!is.na(match(colnames(peakabun), batch_effect[,1]))))
				{
					tkdestroy(dlg)
					tkconfigure(tt,cursor="watch")

					DInd <- gsub("_.*", "", colnames(peakabun))
					DRep = gsub("(.*?)_(.*?)", "\\2", colnames(peakabun))
					if(ctrl_label!=""){
						Factor = Factor[match(gsub("_[0-9]", "", colnames(peakabun)), Factor[,1]),2]
						Factor_label = levels(factor(Factor))
						Factor_label = Factor_label[-which(Factor_label==ctrl_label)]
					} else {
						Factor <- Factor[match(gsub("_[0-9]", "", colnames(peakabun)), Factor[,1]),]
						Factor <- subset(Factor,select=-1) 
						Factor_label <- colnames(Factor)
					}
					Covariate = if(ncol(Covariate)==2){
						as.data.frame(matrix(Covariate[match(gsub("_[0-9]", "", colnames(peakabun)), Covariate[,1]),-1], ncol = ncol(Covariate)-1, dimnames = list(NULL, colnames(Covariate)[2])))
					}else{
						as.data.frame(Covariate[match(gsub("_[0-9]", "", colnames(peakabun)), Covariate[,1]),-1])
					}
					
					for(x in 1:ncol(Covariate)){
						if(cova_type[x]=="D"){
							Covariate[,x] = factor(Covariate[,x])
						}
					}
					colnames(Covariate) = substring(colnames(Covariate), 2, nchar(colnames(Covariate)))
					batch_effect = if(ncol(batch_effect)>2){
						as.data.frame(batch_effect[match(colnames(peakabun), batch_effect[,1]),-1])
					}else{
						if(grepl("latent_group_", tclvalue(textbatchinput))){
							as.data.frame(matrix(factor(batch_effect[match(colnames(peakabun), batch_effect[,1]),-1]), ncol = 1, dimnames = list(NULL, colnames(batch_effect)[2])))
						} else {
							as.data.frame(matrix(batch_effect[match(colnames(peakabun), batch_effect[,1]),-1], ncol = 1, dimnames = list(NULL, colnames(batch_effect)[2])))
						}
					}
					
					
					
					  
					  
					
					peakabun = t(peakabun)
					if(ctrl_label!=""){
						Factor_label = levels(factor(Factor))
						Factor_label = Factor_label[-which(Factor_label==ctrl_label)]
					} else {
						Factor_label <- colnames(Factor)
					}
					
					
					
					
					vif_check <- tktoplevel(); if(isIcon) tk2ico.set(vif_check,icon)
					tkwm.title(vif_check, "Analysis of Covariance - Variance Inflation Factor (VIF)")
					tkpack(fr_head <- tkframe(vif_check), side = "top")
					tkpack(tklabel(fr_head, text = paste("Collinearity Check - VIF", sep = ""), font = fontHeading), side = "top")
					tkpack(fr_vif <- tkframe(vif_check), side = "top")
					
					tclarray <- tclArray()
					VIF_table = c(ifelse(ctrl_label!="Quantity","Group_Comparison",""), "Variable", colnames(batch_effect))
					
					for(lab in Factor_label){
						if(ctrl_label!="") case_idx = which(Factor %in% c(lab, ctrl_label))
						else case_idx = 1:nrow(Factor)
						temp_factor = c(ifelse(ctrl_label!="",paste(lab, "vs.", ctrl_label, sep = ""),lab), ifelse(ctrl_label!="",paste(lab, "-", ctrl_label, "(var)", sep = ""),paste(lab, "(var)", sep = "")), 
						sapply(1:ncol(batch_effect), function(x){
							if(ctrl_label!="") temp = vif(lm(peakabun[case_idx,1]~factor(Factor[case_idx]) + batch_effect[case_idx,x]))
							else temp = vif(lm(peakabun[case_idx,1]~Factor[case_idx,lab] + batch_effect[case_idx,x]))
							ifelse(class(temp)=="matrix", round(temp[2,3], 1), round(temp[2], 1))
						}))
						temp_cova = sapply(1:ncol(Covariate), function(x){
							sapply(1:ncol(batch_effect), function(y){
								temp = vif(lm(peakabun[case_idx,1] ~ Covariate[case_idx,x] + batch_effect[case_idx,y]))
								ifelse(class(temp)=="matrix", round(temp[2,3], 1), round(temp[2], 1))
							})
						})

						temp_cova = rbind(rep(ifelse(ctrl_label!="",paste(lab, "vs.", ctrl_label, sep = ""),lab), ncol(Covariate)), colnames(Covariate), temp_cova)
						VIF_table = cbind(VIF_table, temp_factor, temp_cova)
					}

					for (i in (1:nrow(VIF_table)))
						for (j in (1:ncol(VIF_table)))
							tclarray[[i-1,j-1]] <- as.tclObj(VIF_table[i,j], drop=T)
					scr <- tkscrollbar(fr_vif, repeatinterval = 35, orient="horizontal", command = function(...) tkxview(table1,...))
					table1 <- tkwidget(fr_vif,"table",variable=tclarray,rows=nrow(VIF_table),cols=ncol(VIF_table),titlerows=2,titlecols = 1, selectmode="extended",colwidth=20,background="white", xscrollcommand = function(...) tkset(scr, ...))
					tkpack(table1)
					tkpack(scr, fill = "x")
					tkconfigure(table1,selectmode="extended",rowseparator="\"\n\"",colseparator="\"\t\"")
					tkconfigure(table1,resizeborders="none", state = "disable")    
					tkpack(tklabel(fr_vif, text = "", height = 0, font = fontIntro_para), side = "bottom")
					
					tkpack(tklabel(vif_check, text = "", height = 0, font = fontIntro_para), side = "bottom")
					tkpack(fr_next <- tkframe(vif_check), side = "bottom", fill = "x")
					vif_threshold <- tclVar("10")
					tkpack(tklabel(fr_next,text="       The VIF upper bound:    "), side = "left")
					tkpack(tkentry(fr_next, width="7", textvariable = vif_threshold, bg = "white"), side = "left")
					tkpack(tklabel(fr_next, text = "       "), side = "right")
					

					
					tkpack(tkbutton(fr_next, text = "  Next  ", command = function(...){
						tkdestroy(vif_check)
						tkconfigure(tt,cursor="watch")
						vif_threshold = as.numeric(tclvalue(vif_threshold))
						write.table(VIF_table, file = paste(outpath, "/VIF_threshold_", vif_threshold, ".csv", sep = ""), sep = ",", col.names = F, row.names = F, quote = F)
						
						VIF_table = matrix(as.numeric(VIF_table[-c(1,2),-1]), nrow = nrow(VIF_table)-2, ncol = ncol(VIF_table)-1)
						VIF_table1 = matrix(1, nrow = nrow(VIF_table), ncol = ncol(VIF_table))
						VIF_table1[VIF_table > vif_threshold] = 0
						
						VIF_table = VIF_table1
						
						if(ctrl_label!="")
							ancova_cov = as.vector(sapply(Factor_label, function(x) paste(colnames(Covariate), " (", x, " vs. ", ctrl_label, ")", sep = "")))
						else
							ancova_cov <- as.vector(sapply(Factor_label, function(x) paste(colnames(Covariate), " (", x, ")", sep = "")))
						ancova_var_idx <<- 1:length(ancova_cov)
						ancova_cov_idx <<- NULL

						cova_select <- tktoplevel(); if(isIcon) tk2ico.set(cova_select,icon)
						tkwm.title(cova_select, "Analysis of Covariance - Covariate Selection")
						choose_frame <- tkframe(cova_select)
						var_frame <- tkframe(choose_frame)
						var_scr <- tkscrollbar(var_frame, repeatinterval = 5, command = function(...) tkyview(var_tl,...))
						var_tl <- tklistbox(var_frame, height = 15, selectmode = "extended", yscrollcommand = function(...) tkset(var_scr,...), background = "white")
						tkgrid(tklabel(var_frame, text = "Variable"))
						tkgrid(var_tl, var_scr)
						tkgrid.configure(var_scr, rowspan = 4, sticky = "nsw")
						for(x in ancova_var_idx){
							tkinsert(var_tl, "end", ancova_cov[x])
						}
						cova_frame <- tkframe(choose_frame)
						cova_scr <- tkscrollbar(cova_frame, repeatinterval = 5, command = function(...) tkyview(cova_tl,...))
						cova_tl <- tklistbox(cova_frame, height = 15, selectmode = "extended", yscrollcommand = function(...) tkset(cova_scr, ...), background = "white")
						tkgrid(tklabel(cova_frame, text = "Covariate"))
						tkgrid(cova_tl, cova_scr)
						tkgrid.configure(cova_scr, rowspan = 4, sticky = "nsw")
						in_out <- tkframe(choose_frame)
						tkpack(tkbutton(in_out, text = "  >>  ", command = function(...){
							if(length(as.integer(tkcurselection(var_tl)))!=0){
								varIndex <- as.integer(tkcurselection(var_tl))
								for(x in 1:length(varIndex)){
									tkdelete(var_tl, varIndex[x]-x+1)
								}
								ancova_cov_idx <<- sort(c(ancova_cov_idx, ancova_var_idx[varIndex+1]))
								for(x in varIndex){
									tkinsert(cova_tl, which(ancova_cov_idx==ancova_var_idx[x+1])-1, ancova_cov[ancova_var_idx[x+1]])
								}
								ancova_var_idx <<- ancova_var_idx[-(varIndex+1)]
							}
							tkconfigure(but.covSel,state=ifelse(length(ancova_cov_idx)>0,"normal","disable"))
						}))
						tkpack(tklabel(in_out, text = "     "))
						tkpack(tkbutton(in_out, text = "  <<  ", command = function(...){
							if(length(as.integer(tkcurselection(cova_tl)))!=0){
								covIndex <- as.integer(tkcurselection(cova_tl))
								for(x in 1:length(covIndex)){
									tkdelete(cova_tl, covIndex[x]-x+1)
								}
								ancova_var_idx <<- sort(c(ancova_var_idx, ancova_cov_idx[covIndex+1]))
								for(x in covIndex){
									tkinsert(var_tl, which(ancova_var_idx==ancova_cov_idx[x+1])-1, ancova_cov[ancova_cov_idx[x+1]])
								}
								ancova_cov_idx <<- ancova_cov_idx[-(covIndex+1)]
							}
							tkconfigure(but.covSel,state=ifelse(length(ancova_cov_idx)>0,"normal","disable"))
						}))

						tkgrid(var_frame, in_out, cova_frame, padx = 10)
						tkgrid(tklabel(choose_frame, text = "", height = 0, font = fontIntro_para))
						tkgrid(choose_frame)
						
						
						fr_next <- tkframe(cova_select)
						tkpack(but.covSel <- tkbutton(fr_next, text = "  Next  ", state="disabled", command = function(...){
							tkdestroy(cova_select)
							tkconfigure(tt,cursor="watch")
							if(sum(ancova_cov_idx)){
								ancova_cov[-ancova_cov_idx] = 0
								ancova_cov[ancova_cov_idx] = 1
							}

							if(check) tkmessageBox(title="Analysis of Covariance", message = sprintf("Instead of nested ANCOVA, ANCOVA is applied\nbecause %d samples have only one replicate.",check), icon = "info", type = "ok")

							ancova_cov = matrix(as.numeric(ancova_cov), nrow = length(colnames(Covariate)), dimnames = list(colnames(Covariate), Factor_label))
							nested_ancova <- tktoplevel(); if(isIcon) tk2ico.set(nested_ancova,icon)
							tkwm.title(nested_ancova, "Analysis of Covariance - Model")
							tkpack(fr_head <- tkframe(nested_ancova), side = "top", fill = "x")
							tkpack(tklabel(fr_head, text = "Please input an ANCOVA model."), side = "left")
							tkpack(fr_input <- tkframe(nested_ancova), side = "top", fill = "x")
							for(grp in 1:ncol(ancova_cov)){
								tkpack(fr_input.1 <- tkframe(nested_ancova), side = "top", fill = "x")
								tkpack(fr_input.2 <- tkframe(nested_ancova), side = "top", fill = "x")
								if(all(ancova_cov[,grp]==0)){
									tkpack(tklabel(fr_input.1, text = paste("     Candidate covariates:  ", paste(ifelse(ctrl_label!="","Group , ","Quan , "),ifelse(check,"","Rep , "), paste(rownames(ancova_cov)[which(ancova_cov[,Factor_label[grp]]==1)], collapse = " , "), sep = ""), sep = "")), side = "left")
									if(ctrl_label!="") tkpack(tklabel(fr_input.2, text = paste("     ", Factor_label[grp], " vs. ", ctrl_label, ":   Group + Batch Effect + ", sep = "")), tkentry(fr_input.2, width = "80", textvariable = tclVar("Group:Rep"), bg = "white", state = "disable"), tklabel(fr_input.2, text = "     "), side = "left")
									else tkpack(tklabel(fr_input.2, text = paste("     ", Factor_label[grp], ":   Quan + Batch Effect + ", sep = "")), tkentry(fr_input.2, width = "80", textvariable = tclVar("Quan:Rep"), bg = "white", state = "disable"), tklabel(fr_input.2, text = "     "), side = "left")
								}else{
									tkpack(tklabel(fr_input.1, text = paste("Candidate covariates:  ", paste(ifelse(ctrl_label!="","Group , ","Quan , "),ifelse(check,"","Rep , "), paste(rownames(ancova_cov)[which(ancova_cov[,Factor_label[grp]]==1)], collapse = " , "), sep = ""), sep = ""), wraplength = 750, justify = "left", padx = 17), side = "left")
									assign(paste("model_", Factor_label[grp], sep = ""), tclVar(paste(c(if(check) NULL else ifelse(ctrl_label!="","Group:Rep","Quan:Rep"), rownames(ancova_cov)[which(ancova_cov[,Factor_label[grp]]==1)]), collapse = " + ")))
									if(ctrl_label!="") tkpack(tklabel(fr_input.2, text = paste("     ", Factor_label[grp], " vs. ", ctrl_label, ":   Group + Batch Effect + ", sep = "")), tkentry(fr_input.2, width = "80", textvariable = get(paste("model_", Factor_label[grp], sep = "")), bg = "white"), tklabel(fr_input.2, text = "     "), side = "left")
									else tkpack(tklabel(fr_input.2, text = paste("     ", Factor_label[grp], ":   Quan + Batch Effect + ", sep = "")), tkentry(fr_input.2, width = "80", textvariable = get(paste("model_", Factor_label[grp], sep = "")), bg = "white"), tklabel(fr_input.2, text = "     "), side = "left")
								}
							}
							tkpack(tklabel(nested_ancova, text = "", height = 0, font = fontIntro_para), side = "top")
							tkpack(fr_note <- tkframe(nested_ancova), side = "top", fill = "x")
							tkpack(tklabel(fr_note, text = "", height = 0, font = fontIntro_para), side = "bottom")
							tkpack(tklabel(fr_note, text = "     Note: Notation ''A:B'' indicates that variable B nested in variable A.", justify = "left"), side = "left")
							tkpack(tklabel(fr_note, text = "     "), side = "right")
							tkpack(tkbutton(fr_note, text = "  Next  ", command = function(...){
								model = rep(NA, ncol(ancova_cov))
								for(grp in 1:ncol(ancova_cov)){
									if(all(ancova_cov[,grp]==0)) next
									model[grp] = gsub("\n", "", tclvalue(get(paste("model_", Factor_label[grp], sep = ""))))
								}
								tkdestroy(nested_ancova)
								tkconfigure(tt,cursor="watch")
								op <- options(contrasts = c("contr.helmert", "contr.poly"))
								pb <- tkProgressBar("ANCOVA - Please wait for ANCOVA processing", "0% done", 0, 100, 0, width = 500)
								cat("ANCOVA - Please wait for ANCOVA processing 0 ")
								if(tclvalue(permute.val)=="1" & ncpu>1) cl = makeCluster(ncpu, type = "SOCK")
								for(fa in 1:length(Factor_label)){
									fac = Factor_label[fa]
									check_model = any(ancova_cov[,fac]==1)
									model_choose = model[fa]
									if(ctrl_label!=""){
										group_idx = which(Factor %in% c(fac, ctrl_label))
										Factor_temp = factor(Factor[group_idx])
									} else {
										group_idx <- 1:nrow(Factor)
										Factor_temp <- Factor[,fa]
									}
									Rep = factor(DRep[group_idx])
									Ind = factor(DInd[group_idx],level=unique(DInd[group_idx]))
									
									cov_idx = which(ancova_cov[,fac]==1)
									
									if(check_model){
										be_choose = apply(matrix(VIF_table[,c(1,cov_idx+1)], nrow = nrow(VIF_table)), 1, function(x) all(as.logical(as.numeric(x))))
										
										
										temp_all = data.frame(Covariate[group_idx, cov_idx], batch_effect[group_idx, which(be_choose)])
										colnames(temp_all) = c(colnames(Covariate)[cov_idx], colnames(batch_effect)[which(be_choose)])

										options(contrasts = c("contr.helmert", "contr.poly"))
										anova <- apply(peakabun[group_idx,],2,function(x) {
											nLevel <- apply(cbind(Factor_temp,temp_all)[!is.na(x),],2,function(y) length(unique(y)))
											if(any(nLevel<2)) NA else aov(as.formula(paste("x ~ Factor_temp + ", gsub("Group|Quan", "Factor_temp", model[fa]), "+", paste(colnames(temp_all)[(length(cov_idx)+1):ncol(temp_all)], collapse = "+"))), data = temp_all)
										})
										pvalue <- do.call("rbind",sapply(anova,function(x) if(all(is.na(x))) NA else if(summary(x)[[1]]["Residuals","Sum Sq"]<2.2e-16) NA else Anova(x,type="III",singular.ok=T)[,"Pr(>F)"],simplify=F))
										colnames(pvalue) <- gsub(" ","",gsub("_temp","",gsub("Factor_temp", sprintf("%s%s%s",fac,ifelse(ctrl_label!="","_",""),ctrl_label), c("(Intercept)",rownames(summary(anova[[which(!apply(pvalue,1,function(x) all(is.na(x))))[1]]])[[1]])))))
										options(contrasts = c("contr.SAS", "contr.poly"))
										if(ctrl_label!="") contrasts(Factor_temp) = contrasts(Factor_temp)[c(which(levels(Factor_temp)!=ctrl_label), which(levels(Factor_temp)==ctrl_label)),]
										beta_temp <- apply(peakabun[group_idx,],2,function(x) {
											nLevel <- apply(cbind(Factor_temp,temp_all)[!is.na(x),],2,function(y) length(unique(y)))
											if(any(nLevel<2)) NA else lm(as.formula(paste("x ~ Factor_temp + ", gsub("Group|Quan", "Factor_temp", model[fa]), "+", paste(colnames(temp_all)[(length(cov_idx)+1):ncol(temp_all)], collapse = "+"))), data = temp_all)$coefficients
										})
										beta_temp = switch(class(beta_temp), "matrix" = t(beta_temp), "list"=do.call(rbind, beta_temp))
										resi_temp = apply(peakabun[group_idx,], 2, function(x){
											nLevel <- apply(cbind(Factor_temp,temp_all)[!is.na(x),],2,function(y) length(unique(y)))
											temp = rep(NA, length(x))

											if(all(nLevel>=2)) temp[!is.na(x)&rowSums(is.na(temp_all))==0] = lm(as.formula(paste("x ~ ", gsub("Group|Quan", "Factor_temp", model[fa]), "+", paste(colnames(temp_all)[(length(cov_idx)+1):ncol(temp_all)], collapse = "+"))), data = temp_all)$residuals
											return(temp)
										})
										resi_temp = switch(class(resi_temp), "matrix" = t(resi_temp), "list"=do.call(rbind, resi_temp))
									}else{
										be_choose = as.logical(as.numeric(VIF_table[,1]))
										temp_all = as.data.frame(batch_effect[group_idx, which(be_choose)])
										colnames(temp_all) = colnames(batch_effect)[which(be_choose)]

										options(contrasts = c("contr.helmert", "contr.poly"))
										anova <- apply(peakabun[group_idx,],2,function(x) {
											nLevel <- apply(cbind(Factor_temp,temp_all)[!is.na(x),],2,function(y) length(unique(y)))
											if(any(nLevel<2)) NA else aov(as.formula(paste(ifelse(check,"x ~ Factor_temp + ","x ~ Factor_temp + Factor_temp:Rep + "), paste(colnames(temp_all)[(length(cov_idx)+1):ncol(temp_all)], collapse = "+"))), data = temp_all)
										})
										pvalue <- do.call("rbind",sapply(anova,function(x) if(all(is.na(x))) NA else if(summary(x)[[1]]["Residuals","Sum Sq"]<2.2e-16) NA else Anova(x,type="III",singular.ok=T)[,"Pr(>F)"],simplify=F))
										colnames(pvalue) <- gsub(" ","",gsub("_temp","",gsub("Factor_temp", sprintf("%s%s%s",fac,ifelse(ctrl_label!="","_",""),ctrl_label), c("(Intercept)",rownames(summary(anova[[which(!apply(pvalue,1,function(x) all(is.na(x))))[1]]])[[1]])))))
										options(contrasts = c("contr.SAS", "contr.poly"))
										if(ctrl_label!="") contrasts(Factor_temp) = contrasts(Factor_temp)[c(which(levels(Factor_temp)!=ctrl_label), which(levels(Factor_temp)==ctrl_label)),]
										beta_temp <- apply(peakabun[group_idx,],2,function(x) {
											nLevel <- length(unique(Factor_temp[!is.na(x)]))
											if(any(nLevel<2)) NA else lm(as.formula(paste(ifelse(check,"x ~ Factor_temp + ","x ~ Factor_temp + Factor_temp:Rep + "), paste(colnames(temp_all)[(length(cov_idx)+1):ncol(temp_all)], collapse = "+"))), data = temp_all)$coefficients
										})
										beta_temp = switch(class(beta_temp), "matrix" = t(beta_temp), "list"=do.call(rbind, beta_temp))
										resi_temp = apply(peakabun[group_idx,], 2, function(x){
											nLevel <- apply(cbind(Factor_temp,temp_all)[!is.na(x),],2,function(y) length(unique(y)))
											temp = rep(NA, length(x))
											if(all(nLevel>=2)) temp[!is.na(x)&rowSums(is.na(temp_all))==0] = lm(as.formula(paste(ifelse(check,"x ~ ","x ~ Factor_temp:Rep + "), paste(colnames(temp_all)[(length(cov_idx)+1):ncol(temp_all)], collapse = "+"))), data = temp_all)$residuals
											return(temp)
										})
										resi_temp = switch(class(resi_temp), "matrix" = t(resi_temp), "list"=do.call(rbind, resi_temp))
									}

									pvalue = pvalue[,-ncol(pvalue)]
									write.table(cbind(peak_index, pvalue), file = sprintf("%s/p-value_%s%s%s_details.csv",outpath,fac,ifelse(ctrl_label!=""," vs. ",""),ctrl_label), col.names = T, row.names = F, quote = 1, sep = ",")
									pvalue = as.data.frame(pvalue[,2])
									beta_temp = as.data.frame(beta_temp[,2])
									resi_temp = as.data.frame(resi_temp); colnames(resi_temp) <- rownames(peakabun[group_idx,]) 
									colnames(pvalue) = sprintf("%s%s%s",fac,ifelse(ctrl_label!="","_",""),ctrl_label)
									pfdr = apply(pvalue, 2, function(x) p.adjust(x, method = "BH"))
									colnames(pvalue) = paste(colnames(pvalue), "(pval)", sep = "")
									colnames(pfdr) = paste(colnames(pfdr), "(pfdr)", sep = "")
									colnames(beta_temp) = paste(fac, "_", ctrl_label, "(coef)", sep = "")
									
									write.table(cbind(peak_index, pfdr), file = sprintf("%s/pFDR_%s%s%s.csv",outpath,fac,ifelse(ctrl_label!=""," vs. ",""),ctrl_label), col.names = T, row.names = F, quote = 1, sep = ",")
									write.table(cbind(peak_index, beta_temp), file = sprintf("%s/coef_%s%s%s.csv",outpath,fac,ifelse(ctrl_label!=""," vs. ",""),ctrl_label), col.names = T, row.names = F, quote = 1, sep = ",")
									write.table(cbind(peak_index, resi_temp), file = sprintf("%s/resi_%s%s%s_forPLS.csv",outpath,fac,ifelse(ctrl_label!=""," vs. ",""),ctrl_label), col.names = T, row.names = F, quote = 1, sep = ",")
									if(tclvalue(permute.val)!="1") write.table(cbind(peak_index, resi_temp)[pfdr<.05,], file = sprintf("%s/resi_pFDRsig_%s%s%s_forPLS.csv",outpath,fac,ifelse(ctrl_label!=""," vs. ",""),ctrl_label), col.names = T, row.names = F, quote = 1, sep = ",")

									if(tclvalue(permute.val)=="1"){
										if(sum(pvalue[,1]<as.numeric(tclvalue(permutecut.val)),na.rm=T)>0){
											peak_sig = which(pvalue[,1]<as.numeric(tclvalue(permutecut.val)))
											pvalue_sig = pvalue[peak_sig, 1]
											peakabun_sig = as.data.frame(peakabun[group_idx,peak_sig])
											if(ncpu==1){
												pvalue_permute = sapply(1:as.numeric(tclvalue(permute.times)), function(perm){
													Factor_temp_perm = sample(Factor_temp)
													
													options(contrasts = c("contr.helmert", "contr.poly"))
													if(check_model){
														temp_all2 <- temp_all
														
															
														
															
															

															
															
															
															
														
														anova_temp <- apply(peakabun_sig,2,function(x) {
															nLevel <- length(unique(Factor_temp_perm))
															if(any(nLevel<2)) NA else aov(as.formula(paste("x ~ Factor_temp_perm + ", gsub("Group|Quan", "Factor_temp_perm", model_choose), "+", paste(colnames(temp_all2)[(length(cov_idx)+1):ncol(temp_all2)], collapse = "+"))), data = temp_all2)
														})
														pvalue_temp <- do.call("cbind",sapply(anova_temp,function(x) if(all(is.na(x))) NA else if(summary(x)[[1]]["Residuals","Sum Sq"]<2.2e-16) NA else Anova(x,type="III",singular.ok=T)[,"Pr(>F)"],simplify=F))
														
													}else{
														
															
														
															
															
															
														
														anova_temp <- apply(peakabun_sig,2,function(x) {
															nLevel <- length(unique(Factor_temp_perm))
															if(any(nLevel<2)) NA else aov(as.formula(paste(ifelse(check,"x ~ Factor_temp_perm + ","x ~ Factor_temp_perm + Factor_temp_perm:Rep + "), paste(colnames(temp_all)[(length(cov_idx)+1):ncol(temp_all)], collapse = "+"))), data = temp_all)
															
														})
														pvalue_temp <- do.call("cbind",sapply(anova_temp,function(x) if(all(is.na(x))) NA else if(summary(x)[[1]]["Residuals","Sum Sq"]<2.2e-16) NA else Anova(x,type="III",singular.ok=T)[,"Pr(>F)"],simplify=F))
														
														
													}
													info <- sprintf("%d%%", round(100*fa/length(Factor_label)*(perm/as.numeric(tclvalue(permute.times)))))
													setTkProgressBar(pb, value = round(100*fa/length(Factor_label)*(perm/as.numeric(tclvalue(permute.times)))), sprintf("ANCOVA-Please wait for ANCOVA processing (%s done)", info), paste(fac, "-", ctrl_label, info, sep = " "))
													return(pvalue_temp[2,]<pvalue_sig)
												})
											}else if(ncpu>1){
												clusterExport(cl, c("Factor_temp", "Ind", "Rep", "peakabun_sig", "pvalue_sig", "temp_all", "check_model", "cov_idx", "Anova", "model_choose", "check", "osIsWin"), envir=environment())
												
												pvalue_permute = parSapply(cl, 1:as.numeric(tclvalue(permute.times)), function(perm){

													pb <- tkProgressBar(sprintf("ANCOVA-Parallel task (permutation %d)",perm), "0% done", 0, 100, 0, width = 500)
														Factor_temp_perm = sample(Factor_temp)
														
														options(contrasts = c("contr.helmert", "contr.poly"))
														if(check_model){
															temp_all2 <- temp_all
															
																
															
																
																

																
																
																
																
															
															anova_temp <- apply(peakabun_sig,2,function(x) {
																nLevel <- length(unique(Factor_temp_perm))
																if(any(nLevel<2)) NA else aov(as.formula(paste("x ~ Factor_temp_perm + ", gsub("Group|Quan", "Factor_temp_perm", model_choose), "+", paste(colnames(temp_all2)[(length(cov_idx)+1):ncol(temp_all2)], collapse = "+"))), data = temp_all2)
															})
															pvalue_temp <- do.call("cbind",sapply(anova_temp,function(x) if(all(is.na(x))) NA else if(summary(x)[[1]]["Residuals","Sum Sq"]<2.2e-16) NA else Anova(x,type="III",singular.ok=T)[,"Pr(>F)"],simplify=F))
															
														}else{
															
																
															
																
																
																
															
															anova_temp <- apply(peakabun_sig,2,function(x) {
																nLevel <- length(unique(Factor_temp_perm))
																if(any(nLevel<2)) NA else aov(as.formula(paste(ifelse(check,"x ~ Factor_temp_perm + ","x ~ Factor_temp_perm + Factor_temp_perm:Rep + "), paste(colnames(temp_all)[(length(cov_idx)+1):ncol(temp_all)], collapse = "+"))), data = temp_all)
																
															})
															pvalue_temp <- do.call("cbind",sapply(anova_temp,function(x) if(all(is.na(x))) NA else if(summary(x)[[1]]["Residuals","Sum Sq"]<2.2e-16) NA else Anova(x,type="III",singular.ok=T)[,"Pr(>F)"],simplify=F))
															
															
														}
													if(osIsWin) for(i in 1:100) setTkProgressBar(pb, value = i, sprintf("ANCOVA-Parallel task (permutation %d)", perm), "")
										
													return(pvalue_temp[2,]<pvalue_sig)
												})
											}
											empvalue = as.data.frame((rowSums(pvalue_permute)+1)/(as.numeric(tclvalue(permute.times))+1))
											empvalue = cbind(empvalue, p.adjust(c(empvalue[,1], rep(max(c(empvalue[,1], as.numeric(tclvalue(permutecut.val)))), nrow(peak_index)-nrow(empvalue))), method = "BH")[1:nrow(empvalue)])
											colnames(empvalue) = c(sprintf("%s%s%s(epv)",fac,ifelse(ctrl_label!="","_",""),ctrl_label), sprintf("%s%s%s(epFDR)",fac,ifelse(ctrl_label!="","_",""),ctrl_label))
											write.table(cbind(peak_index[peak_sig,], empvalue), file = sprintf("%s/Empirical p-value_%s%s%s.csv",outpath,fac,ifelse(ctrl_label!=""," vs. ",""),ctrl_label), col.names = T, row.names = F, quote = 1, sep = ",")
											write.table(cbind(peak_index, resi_temp)[peak_sig,][empvalue[,2]<.05,], file = sprintf("%s/resi_EmpiricalpFDRsig_%s%s%s_forPLS.csv",outpath,fac,ifelse(ctrl_label!=""," vs. ",""),ctrl_label), col.names = T, row.names = F, quote = 1, sep = ",")
										}else{
											tkmessageBox(title = "Warning", message = "No significant peaks", icon = "warning", type = "ok")
										}
									}
									
									VIF_table = matrix(VIF_table[,-c(1:(ncol(Covariate)+1))], nrow = nrow(VIF_table))
									info <- sprintf("%d%%", round(100*fa/length(Factor_label)))
									setTkProgressBar(pb, value = round(100*fa/length(Factor_label)), sprintf("ANCOVA - Please wait for ANCOVA processing (%s done)", info), paste(fac, "-", ctrl_label, info, sep = " "))
									if(round(100*fa/length(Factor_label))<100){
										cat(round(100*fa/length(Factor_label)), " ", sep = "")
									}else{
										cat(round(100*fa/length(Factor_label)), " \n", sep = "")
									}
									Sys.sleep(0.1)
								}
								if(tclvalue(permute.val)=="1" & ncpu>1) stopCluster(cl)
								options(op)
								ancova_cov = cbind(c("Cov", rownames(ancova_cov)), rbind(colnames(ancova_cov), ancova_cov))
								write.table(ancova_cov, file = paste(outpath, "/covariate_select.csv", sep = ""), sep = ",", col.names = F, row.names = F, quote = F)
								setTkProgressBar(pb, value = 100, "ANCOVA - Please wait for ANCOVA processing (100% done)", "Finished 100%")
								Sys.sleep(1)
								close(pb)
								tkconfigure(tt,cursor="arrow")
								
								toppeak = tktoplevel(); if(isIcon) tk2ico.set(toppeak,icon)
								tkwm.title(toppeak, "Top peaks")
								fr_top = tkframe(toppeak)
								tkgrid(fr_top)
								top.lab <- tklabel(fr_top, text = "    How many top peaks do you want to paint red?  ")
								top.num = tclVar("50")
								top.entry <- tkentry(fr_top, width = 6, textvariable = top.num, bg = "white",validatecommand="string is integer %P", validate = "all")
								
								tkgrid(top.lab, top.entry, tklabel(fr_top, text = paste("(1 ~ ", ncol(peakabun), ")     ", sep = "")))
								tkgrid(tklabel(toppeak, text = "", height = 0, font = fontIntro_para))
								fr_toppeak = tkframe(toppeak)
								tkgrid(fr_toppeak)
								plot.but <- tkbutton(fr_toppeak, text = "  Plot  ", command = function(...){
									
									tkconfigure(tt,cursor="watch")
									top.num = as.numeric(tclvalue(top.num))
									if(!top.num%in%seq.int(1, ncol(peakabun))) top.num = 50
									for(fa in 1:length(Factor_label)){
										fac = Factor_label[fa]
										pfdr = read.table(sprintf("%s/pFDR_%s%s%s.csv",outpath,fac,ifelse(ctrl_label!=""," vs. ",""),ctrl_label), sep = ",", header = T, quote = "\"")
										pvalue = read.table(sprintf("%s/p-value_%s%s%s_details.csv",outpath,fac,ifelse(ctrl_label!=""," vs. ",""),ctrl_label), sep = ",", header = T, quote = "\"")
										beta_temp = read.table(sprintf("%s/coef_%s%s%s.csv",outpath,fac,ifelse(ctrl_label!=""," vs. ",""),ctrl_label), sep = ",", header = T, quote = "\"")
										
									

										png(filename = sprintf("%s/Volcanoplot_%s%s%s_pFDR_%d.png",outpath,fac,ifelse(ctrl_label!=""," vs. ",""),ctrl_label,top.num), width = 1440, height = 800)
										i=4
										par(mar = c(5,4,1,3))
										pch = rep(21, nrow(pfdr))
										pch[which(-log10(pfdr[,i])>=-log10(0.05))] = 24
										bg = rep("white", nrow(pfdr))
										bg[which(-log10(pfdr[,i])>=-log10(0.05))] = "green3"
										bg[order(-log10(pfdr[,i]), decreasing = T)[1:top.num]] = "red"
										plot(beta_temp[,4], -log10(pfdr[,i]), xlab = "", ylab = "", cex = 1, pch = pch, bg = bg)
										text(beta_temp[order(-log10(pfdr[,i]), decreasing = T)[1:top.num],4], -log10(pfdr[order(-log10(pfdr[,i]), decreasing = T)[1:top.num],i]), label = peak_index[order(-log10(pfdr[,i]), decreasing = T)[1:top.num],1], cex = 1, pos = 1, adj = 0.2)

										abline(h = -log10(0.05), col = "green3")
										par(xpd = NA)
										text(par("usr")[1], -log10(0.05), col = "green3", label = substitute(paste(plain(-log[10]), "(0.05) (", test, " significant peaks)"), list(test = sum(pfdr[,i]<0.05))), font = 4, cex = 1.5, adj = c(0,-0.1))
										par(xpd = F)
										mtext("Fold Change (coefficient)", line = 2, side = 1, col="black", padj = 1, cex = 1.5)
										mtext(expression(-log[10](pFDR)), line = 3.5, side = 2, col="black", padj = 1, cex = 1.5)
										legend(x = "topleft", legend = c("p-FDR<0.05", paste("top", top.num, sep = "")), pch = 24, bty = "o", box.col = "black", cex = 1, pt.bg = c("green3", "red"), bg = rgb(1, 1, 1, alpha = 0.7))
										dev.off()
										
										png(filename = sprintf("%s/Volcanoplot_%s%s%s_pv_%d.png",outpath,fac,ifelse(ctrl_label!=""," vs. ",""),ctrl_label,top.num), width = 1440, height = 800)
										i=5
										par(mar = c(5,4,1,3))
										pch = rep(21, nrow(pvalue))
										pch[which(-log10(pvalue[,i])>=-log10(0.05))] = 24
										bg = rep("white", nrow(pvalue))
										bg[which(-log10(pvalue[,i])>=-log10(0.05))] = "green3"
										bg[order(-log10(pvalue[,i]), decreasing = T)[1:top.num]] = "red"
										plot(beta_temp[,4], -log10(pvalue[,i]), cex = 1, pch = pch, bg = bg, xlab = "", ylab = "")
										text(beta_temp[order(-log10(pvalue[,i]), decreasing = T)[1:top.num],4], -log10(pvalue[order(-log10(pvalue[,i]), decreasing = T)[1:top.num],i]), label = peak_index[order(-log10(pvalue[,i]), decreasing = T)[1:top.num],1], cex = 1, pos = 1, adj = 0.2)
										abline(h = -log10(0.05), col = "green3")
										par(xpd = NA)
										text(par("usr")[1], -log10(0.05), col = "green3", label = substitute(paste(plain(-log[10]), "(0.05) (", test, " significant peaks)"), list(test = sum(pvalue[,i]<0.05))), font = 4, cex = 1.5, adj = c(0,-0.1))
										par(xpd = F)
										mtext("Fold Change (coefficient)", line = 2, side = 1, col="black", padj = 1, cex = 1.5)
										mtext(expression(-log[10](pv)), line = 3.5, side = 2, col="black", padj = 1, cex = 1.5)
										legend(x = "topleft", legend = c("p-value<0.05", paste("top", top.num, sep = "")), pch = 24, bty = "o", box.col = "black", cex = 1, pt.bg = c("green3", "red"), bg = rgb(1, 1, 1, alpha = 0.7))
										dev.off()

										
										if(tclvalue(permute.val)=="1" & file.exists(sprintf("%s/Empirical p-value_%s%s%s.csv",outpath,fac,ifelse(ctrl_label!=""," vs. ",""),ctrl_label))){
											beta_temp = beta_temp[which(pvalue[,5]<as.numeric(tclvalue(permutecut.val))),]
											peak_index_temp = peak_index[which(pvalue[,5]<as.numeric(tclvalue(permutecut.val))),]
											pvalue = read.table(sprintf("%s/Empirical p-value_%s%s%s.csv",outpath,fac,ifelse(ctrl_label!=""," vs. ",""),ctrl_label), sep = ",", header = T, quote = "\"")
											pfdr = pvalue
											
										

											png(filename = sprintf("%s/Volcanoplot_%s%s%s_epFDR_%d.png",outpath,fac,ifelse(ctrl_label!=""," vs. ",""),ctrl_label,top.num), width = 1440, height = 800)
											i=5
											par(mar = c(5,4,1,3))
											pch = rep(21, nrow(pfdr))
											pch[which(-log10(pfdr[,i])>=-log10(0.05))] = 24
											bg = rep("white", nrow(pfdr))
											bg[which(-log10(pfdr[,i])>=-log10(0.05))] = "green3"
											bg[order(-log10(pfdr[,i]), decreasing = T)[1:min(top.num, nrow(pfdr))]] = "red"
											plot(beta_temp[,4], -log10(pfdr[,i]), xlab = "", ylab = "", cex = 1, pch = pch, bg = bg)
											text(beta_temp[order(-log10(pfdr[,i]), decreasing = T)[1:min(top.num, nrow(pfdr))],4], -log10(pfdr[order(-log10(pfdr[,i]), decreasing = T)[1:min(top.num, nrow(pfdr))],i]), label = peak_index_temp[order(-log10(pfdr[,i]), decreasing = T)[1:min(top.num, nrow(pfdr))],1], cex = 1, pos = 1, adj = 0.2)

											abline(h = -log10(0.05), col = "green3")
											par(xpd = NA)
											text(par("usr")[1], -log10(0.05), col = "green3", label = substitute(paste(plain(-log[10]), "(0.05) (", test, " significant peaks)"), list(test = sum(pfdr[,i]<0.05))), font = 4, cex = 1.5, adj = c(0,-0.1))
											par(xpd = F)
											mtext("Fold Change (coefficient)", line = 2, side = 1, col="black", padj = 1, cex = 1.5)
											mtext(expression(-log[10](epFDR)), line = 3.5, side = 2, col="black", padj = 1, cex = 1.5)
											legend(x = "topleft", legend = c("Empirical p-FDR<0.05", paste("top", min(top.num, nrow(pfdr)), sep = "")), pch = 24, bty = "o", box.col = "black", cex = 1, pt.bg = c("green3", "red"), bg = rgb(1, 1, 1, alpha = 0.7))
											dev.off()
											
											png(filename = sprintf("%s/Volcanoplot_%s%s%s_epv_%d.png",outpath,fac,ifelse(ctrl_label!=""," vs. ",""),ctrl_label,top.num), width = 1440, height = 800)
											i=4
											par(mar = c(5,4,1,3))
											pch = rep(21, nrow(pvalue))
											pch[which(-log10(pvalue[,i])>=-log10(0.05))] = 24
											bg = rep("white", nrow(pvalue))
											bg[which(-log10(pvalue[,i])>=-log10(0.05))] = "green3"
											bg[order(-log10(pvalue[,i]), decreasing = T)[1:min(top.num, nrow(pvalue))]] = "red"
											plot(beta_temp[,4], -log10(pvalue[,i]), cex = 1, pch = pch, bg = bg, xlab = "", ylab = "")
											text(beta_temp[order(-log10(pvalue[,i]), decreasing = T)[1:min(top.num, nrow(pvalue))],4], -log10(pvalue[order(-log10(pvalue[,i]), decreasing = T)[1:min(top.num, nrow(pvalue))],i]), label = peak_index_temp[order(-log10(pvalue[,i]), decreasing = T)[1:min(top.num, nrow(pvalue))],1], cex = 1, pos = 1, adj = 0.2)
											abline(h = -log10(0.05), col = "green3")
											par(xpd = NA)
											text(par("usr")[1], -log10(0.05), col = "green3", label = substitute(paste(plain(-log[10]), "(0.05) (", test, " significant peaks)"), list(test = sum(pvalue[,i]<0.05))), font = 4, cex = 1.5, adj = c(0,-0.1))
											par(xpd = F)
											mtext("Fold Change (coefficient)", line = 2, side = 1, col="black", padj = 1, cex = 1.5)
											mtext(expression(-log[10](epv)), line = 3.5, side = 2, col="black", padj = 1, cex = 1.5)
											legend(x = "topleft", legend = c("Empirical p-value<0.05", paste("top", min(top.num, nrow(pvalue)), sep = "")), pch = 24, bty = "o", box.col = "black", cex = 1, pt.bg = c("green3", "red"), bg = rgb(1, 1, 1, alpha = 0.7))
											dev.off()

										}
									}
									tkmessageBox(title = "Analysis of Covariance", message = "Graphs are saved.", icon = "info", type = "ok")
									cat("ANCOVA - Graphs are saved.\n", sep = "")
									tkconfigure(tt,cursor="arrow")
								})
								close.but <- tkbutton(fr_toppeak, text = "  Close  ", command = function(...){
									tkdestroy(toppeak)
									tkmessageBox(title="Analysis of Covariance", message = "ANCOVA is done.", icon = "info", type = "ok")
									cat("ANCOVA - ANCOVA is done.\n", sep = "")
								})
								tkgrid(plot.but, tklabel(fr_toppeak, text = "                              "), close.but)
								tkgrid(tklabel(toppeak, text = "", height = 0, font = fontIntro_para))
							}), side = "right", fill = "y")
							tkconfigure(tt,cursor="arrow")
						}))
						tkgrid(fr_next, sticky = "s", padx = 10)
						tkgrid(tklabel(cova_select, text = "", height = 0, font = fontIntro_para))
						tkmessageBox(title = "ANCOVA - Covariates Select", message = "Please select covariates in ANCOVA.", icon = "info", type = "ok")
						cat("ANCOVA - Please select covariates in ANCOVA.\n")
						tkconfigure(tt,cursor="arrow")
					}), side = "right")
					tkconfigure(tt,cursor="arrow")
					tkmessageBox(title="Analysis of Covariance", message = "Please input the VIF upper bound for claiming a collinearity between batch effect and covariate(s).", icon = "info", type = "ok")
					cat("ANCOVA - Please input the VIF upper bound for claiming a collinearity between batch effect and covariate(s).\n", sep = "")
					
					
					
					
					

				}else{
					tkmessageBox(title = "Error", message = paste("Sample ID not match.\nPlease input the adequate ", paste(c("factor", "covariates", "batch effect")[which(!c(all(!is.na(match(gsub("_[0-9]", "", colnames(peakabun)), Factor[,1]))), all(!is.na(match(gsub("_[0-9]", "", colnames(peakabun)), Covariate[,1]))), all(!is.na(match(colnames(peakabun), batch_effect[,1])))))], collapse = " and "), " file.", sep = ""), icon = "error", type = "ok")
					cat("ANCOVA(Error) - ", paste("Sample ID not match. Please input the adequate ", paste(c("factor", "covariates", "batch effect")[which(!c(all(!is.na(match(gsub("_[0-9]", "", colnames(peakabun)), Factor[,1]))), all(!is.na(match(gsub("_[0-9]", "", colnames(peakabun)), Covariate[,1]))), all(!is.na(match(colnames(peakabun), batch_effect[,1])))))], collapse = " and "), " file", sep = ""), ".\n")
					tkfocus(dlg)
				}
			}
		}else{
			tkmessageBox(title = "Error", message = paste(paste(c("Factor", "Covariate", "Batch effect")[which(c(!file.exists(tclvalue(textFactorinput)), !any(tclvalue(textcova)=="", all(!(tclvalue(textcova)==""), file.exists(tclvalue(textcova)))), !any(tclvalue(textbatchinput)=="", all(!(tclvalue(textbatchinput)==""), file.exists(tclvalue(textbatchinput))))))], collapse = ", "), "file is not found.\nPlease input the correct file path.", sep = ""), icon = "error", type = "ok")
			cat("ANCOVA(Error) - ", paste(paste(c("Factor", "Covariate", "Batch effect")[which(c(!file.exists(tclvalue(textFactorinput)), !any(tclvalue(textcova)=="", all(!(tclvalue(textcova)==""), file.exists(tclvalue(textcova)))), !any(tclvalue(textbatchinput)=="", all(!(tclvalue(textbatchinput)==""), file.exists(tclvalue(textbatchinput))))))], collapse = ", "), "file is not found. Please input the correct file path.", sep = ""), "\n")
			tkfocus(dlg)
		}
	}
	onCancel <- function()
	{
		
		
		tkdestroy(dlg)
		tkfocus(tt)
	}
	fr <- tkframe(dlg)
	OK.but     <-tkbutton(fr,text="     OK     ",command=onOK,state="normal")
	Cancel.but <-tkbutton(fr,text="   Cancel   ",command=onCancel)
	tkgrid(tklabel(fr,text="               "), OK.but,tklabel(fr,text="                                      "), Cancel.but, tklabel(fr,text="               "))
	tkgrid(fr)
	tkgrid(tklabel(dlg, text = "", height = 0, font = fontIntro_para))

	tkfocus(dlg)
	tkbind(dlg, "<Destroy>", function() {tkgrab.release(dlg);tkfocus(tt)})
	
	
	tkwait.window(dlg)

	
}  

PLS <- function(){
	dlg <- tktoplevel(); if(isIcon) tk2ico.set(dlg,icon)
	tkwm.title(dlg, "Partial Least Squares (PLS/PLS-DA)")

	fr_input <- tkframe(dlg)
	tkgrid(tklabel(fr_input, text = "", height = 0, font = fontIntro_para))
	tkgrid(tklabel(fr_input, text = paste("Input file: ", tclvalue(textAbuninput), "    ", sep = ""), padx = 11, wraplength = 600, justify = "left"), sticky = "w")
	tkgrid(tklabel(fr_input, text = paste("Output folder: ", tclvalue(textoutput), "    ", sep = ""), padx = 11, wraplength = 600, justify = "left"), sticky = "w")
	
	fr_input.1 <- tkframe(fr_input)
	tkgrid(fr_input.1, sticky = "w")

	if(!exists("textFactorinput")) textFactorinput <<- tclVar("")
	textFactorWidget <- tkentry(fr_input.1,width="60", textvariable = textFactorinput, bg = "white")
	box.factor <- tkbutton(fr_input.1, text = "...",  command = function() tclvalue(textFactorinput) <- tkgetOpenFile(initialfile = as.character(tclvalue(textFactorinput)), filetypes = "{{Text Files} {.txt .csv}}"))
	fac_label <- tklabel(fr_input.1,text="   Response file (required): ")
	tk2tip(fac_label, "variable of interest")
	tkgrid(fac_label, textFactorWidget, box.factor, tklabel(fr_input.1,text="    "), sticky = "w")
	
	tkgrid(fr_input.1, sticky = "w")
	labelEFactor <- tklabel(fr_input.1, text = "   Number of extracted factors: ")

	textEFactor <- tclVar(2)
	textEFactorWidget <- ttkcombobox(fr_input.1, state = "readonly", values = 2:15, width = 6, textvariable = textEFactor)
	
	
	
	
	tkgrid(labelEFactor, textEFactorWidget, sticky  = "w")
	tkgrid(tklabel(fr_input, text = "", height = 0, font = fontIntro_para))
	tkgrid(fr_input, sticky = "w")

	
	onOK <- function(){
		if(file.exists(tclvalue(textFactorinput))){
			tkconfigure(dlg, cursor = "watch")
			peakabun <- read.table(tclvalue(textAbuninput), header = T, fill = T, sep = ifelse(grepl(".txt", tclvalue(textAbuninput)), "\t", ","), as.is = T, na.string = c("NA", as.character(tclvalue(text_missing))), quote = "\"")  
			resp <- read.table(tclvalue(textFactorinput), header = T, sep = ifelse(grepl(".txt", tclvalue(textFactorinput)), "\t", ","), as.is = T, na.string = c("NA", as.character(tclvalue(text_missing))), quote = "\"")  

			tkconfigure(dlg, cursor = "arrow")


			peak_index <- peakabun[,c(1:3)]
			peakabun <- peakabun[,-c(1:3)]; colnames(peakabun) <- gsub("^\\d{4,8}_(.*?)|(^[A-Z]\\d{4,8})_(.*?)", "\\1", colnames(peakabun))
			resp[,1] <- gsub("^\\d{4,8}_(.*?)|(^[A-Z]\\d{4,8})_(.*?)", "\\1", resp[,1])

			if(all(!is.na(match(gsub("_[0-9]", "", colnames(peakabun)), resp[,1])))){
				tkdestroy(dlg)
				tkconfigure(tt,cursor="watch")

				
				peakabun <- do.call("rbind",by(t(peakabun),gsub("_[0-9]", "", colnames(peakabun)),function(x) apply(x,2,median,na.rm=T)))

				
				peakabun <- peakabun[!apply(peakabun,1,function(x) all(is.na(x))),]

				resp <- resp[match(rownames(peakabun), resp[,1]),]
				resp <- subset(resp,select=-1)
				resp.name <- colnames(resp)

				
				case <- ifelse(ncol(resp)==1,ifelse(toupper(substr(colnames(resp),1,1))=="C",2,1),3) 

				outpath <- paste(tclvalue(textoutput), ifelse(case>1,"/PLS","/PLS-DA"), sep = "")
				dir.create(outpath, showWarnings = F)
				
				switch(as.character(case),	"1"=resp <- as.factor(resp[,1]),
											"2"=resp <- as.numeric(resp[,1]),
											"3"={ resp <- as.matrix(resp)
												  class(resp) <- "numeric" }
				)
				
				if(all(is.na(resp))){
					tkmessageBox(title = "Error", message = "Response file is not correct.", icon = "error", type = "ok")
					cat("PLS/PLS-DA(Error) - Response file is not correct.\n")
					stop("PLS/PLS-DA(Error) - Response file is not correct.\n")
				}

				
				pb <- tkProgressBar("PLS/PLS-DA - Please wait for PLS/PLS-DA processing", "0% done", 0, 100, 0, width = 500)
				cat("PLS/PLS-DA - Please wait for PLS/PLS-DA processing 0 ")

				plsda <- opls(peakabun,resp,predI=as.integer(tclvalue(textEFactor)),orthoI=0,plot=F) 
				
				
				scoreMN <- if(isS4(plsda)) plsda@scoreMN else plsda$scoreMN
				modelDF <- if(isS4(plsda)) plsda@modelDF else plsda$modelDF
				summaryDF <- if(isS4(plsda)) plsda@summaryDF else plsda$summaryDF
				vipVn <- if(isS4(plsda)) plsda@vipVn else plsda$vipVn
				xZeroVarVi <- if(isS4(plsda)) plsda@xZeroVarVi else plsda$xZeroVarVi
				loadingMN <- if(isS4(plsda)) plsda@loadingMN else plsda$loadingMN
				weightMN <- if(isS4(plsda)) plsda@weightMN else plsda$weightMN
				coefficient <- if(isS4(plsda)) plsda@coefficientMN else plsda$coefficients

				
				
					png(sprintf("%s/x_score.png",outpath),1094,1024,res=100)
						plot(plsda,typeVc="x-score",parDevNewL=F)
						
							
							
							
							

							
							

							
							
							
							
								
								
								
								
							
							
								
						
							
							
							
																  
							
							
							
								

							
								
								
								
								
									
									
							
						
					dev.off()
				

				idx <- setdiff(1:nrow(peak_index),xZeroVarVi) 
				
				png(sprintf("%s/VIP%s.png",outpath,ifelse(case!=3,"values","")),1024,1024)
				if(case!=3){ 
					par(mar=c(5,4,4,2)+.1)
					pVn <- apply(peakabun,2,function(x) cor.test(as.integer(as.factor(resp)),x,method=ifelse(case==1,"kendall","pearson"))[["p.value"]])[idx]
					qVn <- qnorm(1 - pVn / 2)
					rmsQuantN <- sqrt(mean(qVn^2))
					par(font=2,font.axis=2,font.lab=2,las=1,mar=c(5.1,4.6,4.1,2.1),lwd=2,pch=16)
					plot(pVn,vipVn,col="#FF000050",pch=16,main="P-value vs. VIP",xlab="p-value",ylab="VIP",xaxs="i",yaxs="i",xlim=c(0,1),cex.axis=1.2,cex.main=1.5,cex.lab=1.4,font.lab=2)
						box()
						curve(qnorm(1 - x / 2) / rmsQuantN,0,1,add=TRUE,col="red",lwd=3)
						abline(h=1,v=0.05,col="blue")

					top10 <- order(vipVn,decreasing=T)[1:10]
					points(pVn[top10],vipVn[top10],col="red",pch=16,xpd=NA)
					text(pVn[top10],vipVn[top10],peak_index[idx,1][top10],col="red",cex=1,font=2,pos=3,xpd=NA)
				} else { 
					top10 <- order(vipVn,decreasing=T)[1:10]
					col.vip <- rep("#FF000050",length(vipVn)); col.vip[top10] <- "red"
					plot(vipVn,col=col.vip,type="h",main="VIP",xlab="",ylab="VIP",xaxs="i",yaxs="i",cex.axis=1.2,cex.main=1.5,cex.lab=1.4,font.lab=2)
						abline(h=1,col="blue")
						points(top10,vipVn[top10],col="red",pch=16,xpd=NA)
						text(top10,vipVn[top10],peak_index[idx,1][top10],col="red",cex=1,font=2,pos=3,xpd=NA)
				}
				dev.off()

				
				confu.tab <- table(resp,fitted(plsda))

				
				write.table(cbind(rownames(scoreMN),scoreMN),sprintf("%s/x_score.csv",outpath),sep=",",quote=F,row.name=F,col.name=c("Sample",colnames(scoreMN)))
				write.table(cbind(peak_index[idx,],loadingMN),sprintf("%s/x_loading.csv",outpath),sep=",",quote=F,row.name=F,col.name=T)
				write.table(cbind(peak_index[idx,],weightMN),sprintf("%s/x_weight.csv",outpath),sep=",",quote=F,row.name=F,col.name=T)
				coeff.name <- if(is.matrix(resp)) resp.name else if (case==1 & length(unique(resp))>2) paste(resp.name,colnames(coefficient),sep="_") else resp.name
				write.table(cbind(peak_index[idx,],coefficient),sprintf("%s/x_coefficient.csv",outpath),sep=",",quote=F,row.name=F,col.name=c(colnames(peak_index),coeff.name))
				write.table(cbind(as.matrix(peak_index[idx,]),vipVn,if(case!=3) pVn[idx]),sprintf("%s/VIP%s.csv",outpath,ifelse(case!=3,"values","")),sep=",",quote=F,row.name=F,col.name=c("Peak_Index","mz","Ret_time.sec","VIP",if(case!=3) "P-value"))

				cat("Variation explanation\n",file=sprintf("%s/summary.csv",outpath))
				write.table(cbind(modelDF,rbind(matrix("",nrow(modelDF)-1,ifelse(case!=3,3,1)),as.matrix(subset(summaryDF,select=c("RMSEE",if(case!=3) c("pR2Y","pQ2")))))),sprintf("%s/summary.csv",outpath),sep=",",quote=F,row.name=F,col.name=T,append=T)
				
					
					
					
					
				

				setTkProgressBar(pb, value = 100, "PLS/PLS-DA - Please wait for PLS/PLS-DA processing (100% done)", "Finished 100%")
				Sys.sleep(1)
				close(pb)
				tkconfigure(tt,cursor="arrow")
				tkmessageBox(title="Partial Least Squares", message = "PLS/PLS-DA is done.", icon = "info", type = "ok")
				cat("PLS/PLS-DA - PLS/PLS-DA is done.\n", sep = "")
			}else{
				tkmessageBox(title = "Error", message = "Sample ID not match.\nPlease input the adequate factors file.", icon = "error", type = "ok")
				cat("PLS/PLS-DA(Error) - Sample ID not match. Please input the adequate factors file.\n")
				tkfocus(dlg)
			}		
		}else{
			tkmessageBox(title = "Error", message = "Response file is not found.\nPlease input the correct file path.", icon = "error", type = "ok")
			cat("PLS/PLS-DA(Error) - ", "Response file is not found. Please input the correct file path.", "\n")
			tkfocus(dlg)
		}
	}
	onCancel <- function()
	{
		tkdestroy(dlg)
		tkfocus(tt)
	}

	fr <- tkframe(dlg)
	OK.but     <-tkbutton(fr,text="     OK     ",command=onOK,state="normal")
	Cancel.but <-tkbutton(fr,text="   Cancel   ",command=onCancel)
	tkgrid(tklabel(fr,text="               "), OK.but,tklabel(fr,text="                                      "), Cancel.but, tklabel(fr,text="               "))
	tkgrid(fr)
	tkgrid(tklabel(dlg, text = "", height = 0, font = fontIntro_para))

	tkfocus(dlg)
	tkbind(dlg, "<Destroy>", function() {tkgrab.release(dlg);tkfocus(tt)})
	
	
	tkwait.window(dlg)
}


xcmsSet = function (files = NULL, snames = NULL, sclass = NULL, phenoData = NULL, 
    profmethod = "bin", profparam = list(), polarity = NULL, 
    lockMassFreq = FALSE, mslevel = NULL, nSlaves = 0, progressCallback = NULL, 
    scanrange = NULL, ...) 
{
    object <- new("xcmsSet")
    xcms.options <- getOption("BioC")$xcms
    xcms.methods <- c(paste("group", xcms.options$group.methods, 
        sep = "."), paste("findPeaks", xcms.options$findPeaks.methods, 
        sep = "."), paste("retcor", xcms.options$retcor.methods, 
        sep = "."), paste("fillPeaks", xcms.options$fillPeaks.methods, 
        sep = "."))
    eval(parse(text = paste("object@progressInfo <- list(", paste(xcms.methods, 
        "=0", sep = "", collapse = ","), ")")))
    if (is.function(progressCallback)) 
        object@progressCallback <- progressCallback
    filepattern <- c("[Cc][Dd][Ff]", "[Nn][Cc]", "([Mm][Zz])?[Xx][Mm][Ll]", 
        "[Mm][Zz][Dd][Aa][Tt][Aa]", "[Mm][Zz][Mm][Ll]")
    filepattern <- paste(paste("\\.", filepattern, "$", sep = ""), 
        collapse = "|")
    if (is.null(files)) 
        files <- getwd()
    info <- file.info(files)
    listed <- list.files(files[info$isdir], pattern = filepattern, 
        recursive = TRUE, full.names = TRUE)
    files <- c(files[!info$isdir], listed)
    files_abs <- file.path(getwd(), files)
    exists <- file.exists(files_abs)
    files[exists] <- files_abs[exists]
    if (lockMassFreq == TRUE) {
        lockMass.files <- grep("02.CDF", files)
        if (length(lockMass.files) > 0) {
            files <- files[-lockMass.files]
        }
    }
    filepaths(object) <- files
    if (length(files) == 0) 
        stop("No NetCDF/mzXML/mzData/mzML files were found.\n")
    fromPaths <- phenoDataFromPaths(files)
    if (is.null(snames)) {
        snames <- rownames(fromPaths)
    }
    else {
        rownames(fromPaths) <- snames
    }
    pdata <- phenoData
    if (is.null(pdata)) {
        pdata <- sclass
        if (is.null(pdata)) 
            pdata <- fromPaths
    }
    phenoData(object) <- pdata
    if (is.null(phenoData)) 
        rownames(phenoData(object)) <- snames
    rtlist <- list(raw = vector("list", length(snames)), corrected = vector("list", 
        length(snames)))
    if ("step" %in% names(profparam)) {
        if ("step" %in% names(list(...)) && profparam$step != 
            list(...)$step) {
            stop("different step values defined in profparam and step arguments")
        }
        profstep <- profparam$step
        profparam <- profparam[names(profparam) != "step"]
    }
    else if ("step" %in% names(list(...))) {
        profstep <- list(...)$step
    }
    else {
        profstep <- 0.1
    }
    if ("method" %in% names(profparam)) {
        if (profparam$method != profmethod) {
            stop("different method values defined in profparam and profmethod arguments")
        }
        profmethod <- profparam$method
        profparam <- profparam[names(profparam) != "method"]
    }
    profinfo(object) <- c(list(method = profmethod, step = profstep), 
        profparam)
    object@polarity <- as.character(polarity)
    includeMSn = FALSE
    includeMSn <- !is.null(mslevel) && mslevel > 1
    xcmsSetArgs <- as.list(match.call())
    if (!is.null(xcmsSetArgs$method)) {
        if (xcmsSetArgs$method == "MS1") {
            includeMSn = TRUE
        }
    }
	
	runParallel <- as.numeric(nSlaves>1)
	parMode <- ifelse(runParallel==1, "SOCK", "")
	snowclust <<- makeCluster(nSlaves, type = "SOCK")
    params <- list(...)
    params$profmethod <- profmethod
    params$profparam <- profparam
    params$includeMSn <- includeMSn
    params$scanrange <- scanrange
    params$mslevel <- mslevel
    params$lockMassFreq <- lockMassFreq
	
    ft <- cbind(file = files, id = 1:length(files))
    argList <- apply(ft, 1, function(x) list(file = x["file"], 
        id = as.numeric(x["id"]), params = params, ratio = as.numeric(x["id"])/nrow(ft)))
    if (parMode == "MPI") {
        res <- xcmsPapply(argList, findPeaksPar)
        mpi.close.Rslaves()
    } else if (parMode == "SOCK") {
		
        res <- xcmsClusterApply(cl = snowclust, x = argList, 
            fun = findPeaksPar, runParallel = runParallel, msgfun = msgfun.featureDetection)
        stopCluster(snowclust)
		rm("snowclust")
    } else {
		
        res <- lapply(argList, findPeaksPar)
		stopCluster(snowclust)
    }
    peaklist <- lapply(res, function(x) x$peaks)
    rtlist$raw <- rtlist$corrected <- lapply(res, function(x) x$scantime)
    if (lockMassFreq) {
        object@dataCorrection[1:length(files)] <- 1
    }
    lapply(1:length(peaklist), function(i) {
        if (is.null(peaklist[[i]])) 
            warning("No peaks found in sample ", snames[i], call. = FALSE)
        else if (nrow(peaklist[[i]]) == 0) 
            warning("No peaks found in sample ", snames[i], call. = FALSE)
        else if (nrow(peaklist[[i]]) == 1) 
            warning("Only 1 peak found in sample ", snames[i], 
                call. = FALSE)
        else if (nrow(peaklist[[i]]) < 10) 
            warning("Only ", nrow(peaklist[[i]]), " peaks found in sample", 
                snames[i], call. = FALSE)
    })
    peaks(object) <- do.call(rbind, peaklist)
    object@rt <- rtlist
    object
}
environment(xcmsSet) <- asNamespace("xcms")
assignInNamespace("xcmsSet", xcmsSet, ns = "xcms", envir = as.environment(asNamespace("xcms")))


xcmsClusterApply_my <- function(cl, x, fun, runParallel, msgfun=NULL, ...) {
    argfun <- function(i) c(list(x[[i]]), list(...))
    n <- length(x)
    checkCluster(cl)
    p <- length(cl)
    if (n > 0 && p > 0) {
        submit <- function(node, job) sendCall(cl[[node]], fun,
                                               argfun(job), tag = job)
        for (i in 1:min(n, p)) {
            if (!is.null(msgfun))
                do.call(msgfun,args=list(x=x,i=i));

            submit(i, i)
			if(runParallel==1){
				info <- sprintf("%d%%", round(100*x[[i]]$ratio))
				
				
				setTkProgressBar(pb, value = round(100*x[[i]]$ratio), title = sprintf("Peak Alignment and Annotation - Please wait for alignment processing (%s done)", info), label = gsub(".*/(.*?).mzXML", "\\1", x[[i]]$file))
				if(round(100*x[[i]]$ratio)<100){
					cat(round(100*x[[i]]$ratio), " ", sep = "")
				}else{
					cat(round(100*x[[i]]$ratio), " \n", sep = "")
				}
			}
        }
        val <- vector("list", n)
        for (i in 1:n) {
            d <- recvOneResult(cl)
            j <- i + min(n, p)
            if (j <= n) {
                if (!is.null(msgfun))
                    do.call(msgfun,args=list(x=x,i=j));
				if(runParallel==1){
					info <- sprintf("%d%%", round(100*x[[j]]$ratio))
					
					
					setTkProgressBar(pb, value = round(100*x[[j]]$ratio), title = sprintf("Peak Alignment and Annotation - Please wait for alignment processing (%s done)", info), label = gsub(".*/(.*?).mzXML", "\\1", x[[j]]$file))
					if(round(100*x[[j]]$ratio)<100){
						cat(round(100*x[[j]]$ratio), " ", sep = "")
					}else{
						cat(round(100*x[[j]]$ratio), " \n", sep = "")
					}
				}
                submit(d$node, j)
            }
            val[d$tag] <- list(d$value)
        }
        checkForRemoteErrors(val)
    }

}
assignInNamespace("xcmsClusterApply", xcmsClusterApply_my, ns="xcms")


findPeaksPar_my <- function(arg) 
{
    require(xcms)
    params <- arg$params
    myID <- arg$id
    if (is.null(params$method)) 
		params$method <- getOption("BioC")$xcms$findPeaks.method
    method <- match.arg(params$method, getOption("BioC")$xcms$findPeaks.methods)
    if (is.na(method)) 
		stop("unknown method : ", method)
    method <- paste("findPeaks", method, sep = ".")
    xRaw <- xcmsRaw(arg$file, profmethod = params$profmethod, 
		profparam = params$profparam, profstep = 0, includeMSn = params$includeMSn, 
		mslevel = params$mslevel, scanrange = params$scanrange)
    if (params$lockMassFreq == TRUE) {
		xRaw <- stitch(xRaw, AutoLockMass(xRaw))
    }
    params["object"] <- xRaw
    params$scanrange <- switch(params$lockMassFreq, T = seq.int(1, length(xRaw@scanindex)-1), F = NULL)
    params["method"] <- params["id"] <- params["profmethod"] <- params["profparam"] <- params["includeMSn"] <- params["lockMassFreq"] <- params["mslevel"] <- NULL
    peaks <- do.call(method, params)
	if(exists("snowclust")){
		info <- sprintf("%d%%", round(100*arg$ratio))
		
		
		setTkProgressBar(pb, value = round(100*arg$ratio), title = sprintf("Peak Alignment and Annotation - Please wait for alignment processing (%s done)", info), label = gsub(".*/(.*?).mzXML", "\\1", arg$file))
		if(round(100*arg$ratio)<100){
			cat(round(100*arg$ratio), " ", sep = "")
		}else{
			cat(round(100*arg$ratio), " \n", sep = "")
		}
	}
    list(scantime = xRaw@scantime, peaks = cbind(peaks, sample = rep.int(myID, nrow(peaks))))
}
assignInNamespace("findPeaksPar", findPeaksPar_my, ns="xcms")



setGeneric("retcor.peakgroups", function(object, ...) standardGeneric("retcor.peakgroups"))
setMethod("retcor.peakgroups", "xcmsSet", function(object, missing = 1, extra = 1,
                                                   smooth = c("loess", "linear"), span = .2,
                                                   family = c("gaussian", "symmetric"),
                                                   plottype = c("none", "deviation", "mdevden"),
                                                   col = NULL, ty = NULL) {

    peakmat <- peaks(object)
    groupmat <- groups(object)
    if (length(groupmat) == 0){
		cat("Alignment (rector.peakgroups error) - No group information found.\n")
		tkmessageBox(title = "Error", message = "No group information found.", icon = "error", type = "ok")
		stop("No group information found")
	}
    samples <- sampnames(object)
    classlabel <- as.vector(unclass(sampclass(object)))
    n <- length(samples)
    corpeaks <- peakmat
    smooth <- match.arg(smooth)
    plottype <- match.arg(plottype)
    family <- match.arg(family)
    if (length(object@rt) == 2)
        rtcor <- object@rt$corrected
    else {
        fnames <- filepaths(object)
        rtcor <- vector("list", length(fnames))
        for (i in seq(along = fnames)) {
            xraw <- xcmsRaw(fnames[i])
            rtcor[[i]] <- xraw@scantime
        }
        object@rt <- list(raw = rtcor, corrected = rtcor)
    }

    nsamp <- rowSums(groupmat[,match("npeaks", colnames(groupmat))+unique(classlabel),drop=FALSE])

    idx <- which(nsamp >= n-missing & groupmat[,"npeaks"] <= nsamp + extra)
    if (length(idx) == 0){
		cat("Alignment (rector.peakgroups error) - No peak groups found for retention time correction.\n")
		tkmessageBox(title = "Error", message = "No peak groups found for retention time correction.", icon = "error", type = "ok")
		stop("No peak groups found for retention time correction")
	}
    idx <- idx[order(groupmat[idx,"rtmed"])]

    rt <- groupval(object, "maxint", "rt")[idx,, drop=FALSE]
    cat("Retention Time Correction Groups:", nrow(rt), "\n")
    rtdev <- rt - apply(rt, 1, median, na.rm = TRUE)

    if (smooth == "loess") {
        mingroups <- min(colSums(!is.na(rt)))
        if (mingroups < 4) {
            smooth <- "linear"
            warning("Too few peak groups, reverting to linear method")
        } else if (mingroups*span < 4) {
            span <- 4/mingroups
            warning("Span too small, resetting to ", round(span, 2))
        }
    }

    rtdevsmo <- vector("list", n)

    
    rtdevrange <- range(rtdev, na.rm = TRUE)
    warn.overcorrect <- FALSE

    for (i in 1:n) {

        pts <- na.omit(data.frame(rt = rt[,i], rtdev = rtdev[,i]))

        if (smooth == "loess") {
            lo <- suppressErrors(loess(rtdev ~ rt, pts, span = span, degree = 1, family = family))

            rtdevsmo[[i]] <- xcms:::na.flatfill(predict(lo, data.frame(rt = rtcor[[i]])))

            rtdevsmo[[i]][abs(rtdevsmo[[i]]) > quantile(abs(rtdevsmo[[i]]), 0.9)*2] <- NA

            if (length(naidx <- which(is.na(rtdevsmo[[i]]))))
                rtdevsmo[[i]][naidx] <- suppressErrors(approx(na.omit(data.frame(rtcor[[i]], rtdevsmo[[i]])),
                                                                xout = rtcor[[i]][naidx], rule = 2)$y)
            while (length(decidx <- which(diff(rtcor[[i]] - rtdevsmo[[i]]) < 0))) {
                d <- diff(rtcor[[i]] - rtdevsmo[[i]])[tail(decidx, 1)]
                rtdevsmo[[i]][tail(decidx, 1)] <- rtdevsmo[[i]][tail(decidx, 1)] - d
                if (abs(d) <= 1e-06)
                    break;
            }

            rtdevsmorange <- range(rtdevsmo[[i]])
            if (any(rtdevsmorange/rtdevrange > 2)) warn.overcorrect <- TRUE
        } else {
            if (nrow(pts) < 2) {
                stop("Not enough ``well behaved'' peak groups even for linear smoothing of retention times")
            }
            fit <- lsfit(pts$rt, pts$rtdev)
            rtdevsmo[[i]] <- rtcor[[i]] * fit$coef[2] + fit$coef[1]
            ptsrange <- range(pts$rt)
            minidx <- rtcor[[i]] < ptsrange[1]
            maxidx <- rtcor[[i]] > ptsrange[2]
            rtdevsmo[[i]][minidx] <- rtdevsmo[[i]][head(which(!minidx), n = 1)]
            rtdevsmo[[i]][maxidx] <- rtdevsmo[[i]][tail(which(!maxidx), n = 1)]
        }
    }

    if (warn.overcorrect) {
        warning(paste("Fitted retention time deviation curves exceed points by more than 2x.",
                      "This is dangerous and the algorithm is probably overcorrecting your data.",
                      "Consider increasing the span parameter or switching to the linear smoothing method.",
                      sep = "\n"))
    }

    if (plottype == "mdevden") {
        split.screen(matrix(c(0, 1, .3, 1, 0, 1, 0, .3), ncol = 4, byrow = TRUE))
        screen(1)
        par(mar = c(0, 4.1, 4.1, 2), xaxt = "n")
    }

    if (plottype %in% c("deviation", "mdevden")) {


        if (missing(col) || is.null(col)) {
            col <- integer(n)
            for (i in 1:max(classlabel))
                col[classlabel == i] <- 1:sum(classlabel == i)
        }
        if (missing(ty) || is.null(ty)) {
            ty <- integer(n)
            for (i in 1:max(col))
                ty[col == i] <- 1:sum(col == i)
        }
        if (length(palette()) < max(col))
            mypal <- rainbow(max(col), end = 0.85)
        else
            mypal <- palette()[1:max(col)]

        rtrange <- range(do.call(c, rtcor))
        devrange <- range(do.call(c, rtdevsmo))

        plot(0, 0, type="n", xlim = rtrange, ylim = devrange, main = "Retention Time Deviation vs. Retention Time", xlab = "Retention Time", ylab = "Retention Time Deviation")
        

        for (i in 1:n) {
            points(data.frame(rt = rt[,i], rtdev = rtdev[,i]), col = mypal[col[i]], pch = ty[i], type="p")
            points(rtcor[[i]], rtdevsmo[[i]], type="l", col = mypal[col[i]], lty = ty[i])
        }
    }

    if (plottype == "mdevden") {

        screen(2)
        par(mar = c(5.1, 4.1, 0, 2), yaxt = "n")
        allden <- density(peakmat[,"rt"], bw = diff(rtrange)/200, from = rtrange[1], to = rtrange[2])[c("x","y")]
        corden <- density(rt, bw = diff(rtrange)/200, from = rtrange[1], to = rtrange[2], na.rm = TRUE)[c("x","y")]
        allden$y <- allden$y / sum(allden$y)
        corden$y <- corden$y / sum(corden$y)
        maxden <- max(allden$y, corden$y)
        plot(c(0,0), xlim = rtrange, ylim = c(0, maxden), type = "n", main = "", xlab = "Retention Time", ylab = "Peak Density")
        points(allden, type = "l", col = 1)
        points(corden, type = "l", col = 2)
        abline(h = 0, col = "grey")
        legend(x = "topright", maxden, c("All", "Correction"), col = 1:2, lty = c(1,1), xjust = 1)
        close.screen(all = TRUE)
    }

    for (i in 1:n) {

        cfun <- stepfun(rtcor[[i]][-1] - diff(rtcor[[i]])/2, rtcor[[i]] - rtdevsmo[[i]])
        rtcor[[i]] <- rtcor[[i]] - rtdevsmo[[i]]

        sidx <- which(corpeaks[,"sample"] == i)
        corpeaks[sidx, c("rt", "rtmin", "rtmax")] <- cfun(corpeaks[sidx, c("rt", "rtmin", "rtmax")])
    }

    object@rt$corrected <- rtcor
    peaks(object) <- corpeaks
    groups(object) <- matrix(nrow = 0, ncol = 0)
    groupidx(object) <- list()
    invisible(object)
})

align_loop <- function(loop, xset_g_rt, outpath, bw, minfrac, minsamp, mzwid, max, missing_n, extra, smooth, span, family, fillpeak, PAformat, IP_method, filter.iso, filter.frg, filter.add, add.pn){
	
	xset_g <- group(xset_g_rt, mzwid = mzwid, minsamp = minsamp, bw = bw, minfrac = minfrac, max = max)
	group_retcor <<- tktoplevel(); if(isIcon) tk2ico.set(group_retcor,icon)
	temp = c("second", "third", "forth", "fifth (final)")
	tkwm.title(group_retcor, paste("Align retention times across samples (", temp[loop-1]," loop)", sep = ""))
	fr_loop <- tkframe(group_retcor)
	png(filename = paste(outpath, "/temp.png", sep = ""), width = 600, height = 450)
	xset_g_rt <- retcor(xset_g, missing = missing_n, extra = extra, smooth = smooth, plottype="m", span = span, family = family)
	dev.off()
	if(exists("imgrvalue")) tkpack.forget(imgrvalue)
	imgrvalue <<- tclVar()
	tkimage.create("photo", imgrvalue, file = paste(outpath, "/temp.png", sep = ""))
	imgAsLabel <- tklabel(fr_loop, image = imgrvalue)
	tkgrid(imgAsLabel, sticky = "w")
	tkgrid(tklabel(fr_loop, text="    ", height = 0, font = fontIntro_para))
	tkgrid(fr_loop, sticky = "w")
	fr_conti <- tkframe(group_retcor)
	Stop.but <- tkbutton(fr_conti,text="   Stop   ",command=function(...){
		tkdestroy(group_retcor)
		tkconfigure(tt,cursor="watch")
		xset_final = group(xset_g_rt, mzwid = mzwid, minsamp = minsamp, bw = bw, minfrac = minfrac, max = max)
		if(fillpeak==1){
			xset_final <- fillPeaks(xset_final, method="chrom")
		}

		
		xset_peakTable = peakTable(xset_final, value = list(c("into", "intf", "maxo", "maxf"), c("into", "intb", "maxo"))[[IP_method]][PAformat])
		xset_peakTable = data.frame(Peak_Index = 1:nrow(xset_peakTable), mz = xset_peakTable[,1], "Ret_Time.sec" = round(xset_peakTable[,4], 4), xset_peakTable[,9:ncol(xset_peakTable)])
		colnames(xset_peakTable)[4:ncol(xset_peakTable)] = gsub("^\\d{4,8}_(.*?)|(^[A-Z]\\d{4,8})_(.*?)", "\\1", colnames(xset_peakTable)[4:ncol(xset_peakTable)])
		write.table(xset_peakTable, file = paste(outpath, "/S", ncol(xset_peakTable)-3, ifelse(fillpeak==1, "_filledPeak", "_nofilledPeak"), list(c("into", "intf", "maxo", "maxf"), c("into", "intb", "maxo"))[[IP_method]][PAformat], "_loop", loop, ".csv", sep = ""), col.names = T, row.names = F, quote = F, sep = ",")

			
			xsa <- xsAnnotate(xset_final)
			
			xsaF <- groupFWHM(xsa)
			
			
			xsaC <- groupCorr(xsaF)
			
			xsaFI <- findIsotopes(xsaC)
			
			xsaFA <- findAdducts(xsaFI, polarity=switch(add.pn,"positive","negative"))
			
			peakList <- getPeaklist(xsaFA)

		
		xset_peakTable <- peakList
		xset_peakTable <- data.frame(Peak_Index = 1:nrow(xset_peakTable), mz = xset_peakTable[,1], "Ret_Time.sec" = round(xset_peakTable[,4], 4), xset_peakTable[,9:ncol(xset_peakTable)])
		colnames(xset_peakTable)[4:ncol(xset_peakTable)] = gsub("^\\d{4,8}_(.*?)|(^[A-Z]\\d{4,8})_(.*?)", "\\1", colnames(xset_peakTable)[4:ncol(xset_peakTable)])
		write.table(xset_peakTable[,c(1:3,(ncol(xset_peakTable)-2):ncol(xset_peakTable))], file = paste(outpath, "/S", ncol(xset_peakTable)-6, ifelse(fillpeak==1, "_filledPeak", "_nofilledPeak"), list(c("into", "intf", "maxo", "maxf"), c("into", "intb", "maxo"))[[IP_method]][PAformat], "_loop", loop, "_annoinfo", ".csv", sep = ""), col.names = T, row.names = F, quote = F, sep = ",")

		if(any(c(filter.iso,filter.frg,filter.add)==1)){
			idx.frg <- idx.iso <- idx.add <- 1:nrow(peakList)
			if(filter.frg==1){
				avgPeakAbun <- rowMeans(peakList[,10:(ncol(peakList)-3)],na.rm=T); names(avgPeakAbun) <- 1:nrow(peakList)
				idx.frg <- as.integer(tapply(avgPeakAbun,peakList$pcgroup,function(x) names(x)[which.max(x)]))
			}
			if(filter.iso==1){
				isotope <- gsub("\\[[0-9]+\\]","",peakList$isotopes)
				idx.iso <- which(isotope%in%c("[M]+",""))
			}
			if(filter.add==1){
				adduct <- strsplit(gsub(" [0-9.]+","",peakList$adduct)," ")
				add.sel <<- add <- sort(unique(unlist(adduct))); add <- add[add!=""]

				
				adduct_sel <<- tktoplevel(); if(isIcon) tk2ico.set(adduct_sel,icon)
				tkwm.title(adduct_sel,"Adduct selection:")
				fr <- tkframe(adduct_sel)
				fr.1 <- tkframe(fr)
					tkgrid.configure(tklabel(fr.1,text="Select adducts to reserve:"),sticky="n")
				fr.2 <- tkframe(fr)
					add.ls <- tclVar(""); tclvalue(add.ls) <- add
					add.scr <- tkscrollbar(fr.2,repeatinterval=10,command=function(...)tkyview(add.lst,...))
					add.lst <- tklistbox(fr.2,width=18,height=10,selectmode="multiple",listvariable=add.ls,yscrollcommand=function(...)tkset(add.scr,...),background="white",cursor="hand2")
				tkgrid(add.lst,add.scr)
				OK.but <- tkbutton(fr,text="     OK     ",command=function(x){
					add.sel <<- add[as.numeric(tkcurselection(add.lst))+1]
					tkdestroy(adduct_sel)
				},state="normal")
				tkgrid.configure(add.scr,sticky="ns")
				tkgrid(tklabel(fr,text=""))
				tkgrid(fr.1)
				tkgrid(fr.2)
				tkgrid(tklabel(fr,text=""))
				tkgrid(OK.but)
				tkgrid(tklabel(fr,text=""))
				tkgrid(tklabel(adduct_sel,text="",width=5),fr,tklabel(adduct_sel,text="",width=5))

				tkwait.window(adduct_sel)

				idx.add <- which(sapply(adduct,function(a) any(a%in%add.sel)) | sapply(adduct,length)==0)
			}
			idx <- sort(intersect(intersect(idx.frg,idx.iso),idx.add))

			
			xset_peakTable <- xset_peakTable[idx,]
			write.table(xset_peakTable[,1:(ncol(xset_peakTable)-3)], file = paste(outpath, "/S", ncol(xset_peakTable)-6, ifelse(fillpeak==1, "_filledPeak", "_nofilledPeak"), list(c("into", "intf", "maxo", "maxf"), c("into", "intb", "maxo"))[[IP_method]][PAformat], "_loop", loop, "_filterAnnoPeak", ".csv", sep = ""), col.names = T, row.names = F, quote = F, sep = ",")
			write.table(xset_peakTable[,c(1:3,(ncol(xset_peakTable)-2):ncol(xset_peakTable))], file = paste(outpath, "/S", ncol(xset_peakTable)-6, ifelse(fillpeak==1, "_filledPeak", "_nofilledPeak"), list(c("into", "intf", "maxo", "maxf"), c("into", "intb", "maxo"))[[IP_method]][PAformat], "_loop", loop, "_filterAnnoPeak_annoinfo", ".csv", sep = ""), col.names = T, row.names = F, quote = F, sep = ",")
		}
		cat("Peak Alignment and Annotation - Peak Annotation is done.\n")

		
		textAbuninput <<- tclVar(paste(outpath, "/S", ncol(xset_peakTable)-6, ifelse(fillpeak==1, "_filledPeak", "_nofilledPeak"), list(c("into", "intf", "maxo", "maxf"), c("into", "intb", "maxo"))[[IP_method]][PAformat], "_loop", loop, ifelse(any(c(filter.iso,filter.frg,filter.add)==1), "_filterAnnoPeak", ""), ".csv", sep = ""))
		unlink(paste(outpath, "/temp.png", sep = ""))
		text_zero <<- tclVar("0")
		text_missing <<- tclVar("NA")
		tkconfigure(tt,cursor="arrow")
		tkmessageBox(title = "Peak Alignment and Annotation", message = "Peak Alignment and Annotation is done.", icon = "info", type = "ok")
		cat("Peak Alignment and Annotation - Peak Alignment and Annotation is done.\n")
		tkentryconfigure(peakMenu, 1, state = "active")
		tkentryconfigure(QCsubMenu, 0, state = "active")
		
		tkentryconfigure(batchsubMenu, 0, state = "active")
		tkentryconfigure(batchsubMenu, 1, state = "active")
		tkentryconfigure(statMenu, 1, state = "active")
		tkentryconfigure(statMenu, 2, state = "active")
		
	})
	if(loop<5){
		Next.but <- tkbutton(fr_conti,text="  Correct again   ",command=function(...){
			tkdestroy(group_retcor)
			tkconfigure(tt,cursor="watch")
			unlink(paste(outpath, "/temp.png", sep = ""))
			align_loop(loop+1, xset_g_rt, outpath, bw, minfrac, minsamp, mzwid, max, missing_n, extra, smooth, span, family, fillpeak, PAformat, IP_method, filter.iso, filter.frg, filter.add, add.pn)
			tkconfigure(tt,cursor="arrow")
		})
		tkgrid(tklabel(fr_conti,text="                    "), Next.but,tklabel(fr_conti,text="                                                "), Stop.but, tklabel(fr_conti,text="                    "))
	}else{
		tkgrid(Stop.but)
	}
	tkgrid(fr_conti)
	tkgrid(tklabel(group_retcor, text="    ", height = 0, font = fontIntro_para))
	tkfocus(group_retcor)
}




Align <- function(ALtitle, afterQC = F, QCpath)  
{
	citation("xcms")
	citation("CAMERA")
	xcmscitation <- tktoplevel(); if(isIcon) tk2ico.set(xcmscitation,icon)
	tkwm.title(xcmscitation, "Citations of XCMS and CAMERA")
	tkpack(tklabel(xcmscitation,text="", height = 0, font = fontIntro_para))
	tkpack(tklabel(xcmscitation, text = "XCMS", font=fontHeading, justify = "left"))
	tkpack(tklabel(xcmscitation, text = "
	Smith CA, Want EJ, O'Maille G, Abagyan R, Siuzdak G. XCMS: processing mass spectrometry data for metabolite
	profiling using nonlinear peak alignment, matching, and identification. Anal Chem. 2006;78(3):779-87.

	Tautenhahn R, Bottcher C, Neumann S. Highly sensitive feature detection for high resolution LC/MS. BMC 
	Bioinformatics. 2008;9:504.

	Benton HP, Want EJ, Ebbels TM. Correction of mass calibration gaps in liquid chromatography-mass spectrometry 
	metabolomics data. Bioinformatics. 2010;26(19):2488-9.
	", justify = "left"))
	tkpack(tklabel(xcmscitation,text="", height = 0, font = fontIntro_para))
	tkpack(tklabel(xcmscitation, text = "CAMERA", font=fontHeading, justify = "left"))
	tkpack(tklabel(xcmscitation, text = "
	Kuhl C, Tautenhahn R, Bottcher C, Larson TR, Neumann S. CAMERA: an integrated strategy for compound spectra   
	extraction and annotation of liquid chromatography/mass spectrometry data sets. Anal Chem. 2012;84(1):283-9.
	", justify = "left"))
	tkpack(tkbutton(xcmscitation, text = "  Close  ", command = function() tkdestroy(xcmscitation)))
	tkpack(tklabel(xcmscitation,text="", height = 0, font = fontIntro_para))
	dlg <- tktoplevel(width = 800); if(isIcon) tk2ico.set(dlg,icon)
	tkwm.title(dlg, ALtitle)
	
	fr_input <- tkframe(dlg)
	tkgrid(tklabel(fr_input,text="", height = 0, font = fontIntro_para))
	if(!exists("textmzdatainput")){
		textmzdatainput <<- tclVar()
		textinputWidget <- tkentry(fr_input,width="70", textvariable = textmzdatainput, bg = "white")
		box.input <- tkbutton(fr_input, text = "...",  command = function(){
			tclvalue(textmzdatainput) <- tkchooseDirectory()
		})
		tkgrid(tklabel(fr_input, text = "   Input folder (mzXML file): "), textinputWidget, box.input, tklabel(fr_input, text = "    "), sticky = "w")
	}else{
		if(tclvalue(textmzdatainput)==""){
			textmzdatainput <<- tclVar()
			textinputWidget <- tkentry(fr_input,width="70", textvariable = textmzdatainput, bg = "white")
			box.input <- tkbutton(fr_input, text = "...",  command = function(){
				tclvalue(textmzdatainput) <- tkchooseDirectory()
			})
			tkgrid(tklabel(fr_input, text = "   Input folder (mzXML file): "), textinputWidget, box.input, tklabel(fr_input, text = "    "), sticky = "w")
		}else{
			tkgrid(tklabel(fr_input, text = paste("Input folder: ", tclvalue(textmzdatainput), "    ", sep = ""), padx = 11, wraplength = 720, justify = "left"), sticky = "w")
		}
	}
	
	if(!exists("textoutput")){
		textoutput <<- tclVar()
		textoutputWidget <- tkentry(fr_input,width="70", textvariable = textoutput, bg = "white")
		box.output <- tkbutton(fr_input, text = "...",  command = function(){
			tclvalue(textoutput) <- tkchooseDirectory()
		})
		tkgrid(tklabel(fr_input, text = "   Output folder: "), textoutputWidget, box.output, tklabel(fr_input, text = "    "), sticky = "w")
	}else{
		if(tclvalue(textoutput)==""){
			textoutput <<- tclVar()
			textoutputWidget <- tkentry(fr_input,width="70", textvariable = textoutput, bg = "white")
			box.output <- tkbutton(fr_input, text = "...",  command = function(){
				tclvalue(textoutput) <- tkchooseDirectory()
			})
			tkgrid(tklabel(fr_input, text = "   Output folder: "), textoutputWidget, box.output, tklabel(fr_input, text = "    "), sticky = "w")
		}else{
			if(tclvalue(textmzdatainput)==""){
				tkgrid(tklabel(fr_input, text = "   Output folder: "), tklabel(fr_input, text = tclvalue(textoutput)), tklabel(fr_input, text = "    "), tklabel(fr_input, text = "    "), sticky = "w")
			}else{
				tkgrid(tklabel(fr_input, text = paste("Output folder: ", tclvalue(textoutput), "    ", sep = ""), padx = 11, wraplength = 720, justify = "left"), sticky = "w")
			}
		}
	}
	
	
	tkgrid(tklabel(fr_input,text="", height = 0, font = fontIntro_para))
	tkpack(fr_input, side = "top", fill = "x")
	
	
	fr_para <- tkframe(dlg)
	tkpack(fr_para, side = "top", fill = "x")
	tkgrid(tklabel(fr_para, text = "   Parameter setting:"))
	
	
	fr_IP <- tkframe(dlg)
	tkpack(fr_IP, side = "top", fill = "x")
	

	
	fr_IP_head = tkframe(fr_IP)
	tkgrid(tklabel(fr_IP_head, text = "         Identify peaks: "))
	tkgrid(fr_IP_head, sticky = "w")
	
	fr_IP_para <- tkframe(fr_IP)
	tkgrid(fr_IP_para, sticky = "w")
	
	labellock <- tklabel(fr_IP_para, text = "                  lockMassFreq: ")
	tk2tip(labellock, "XCMS: Performs correction for Waters LockMass function")
	textlock <- tclVar("F")
	textlockWidget <- tkentry(fr_IP_para, width = 4, textvariable = textlock, bg = "white")
	
	labelpolar <- tklabel(fr_IP_para, text = " , polarity: ")
	tk2tip(labelpolar, "XCMS: filter raw data for positive(+)/negative(-)/both(b) scans")
	textpolar <- tclVar("b")
	textpolarWidget <- tkentry(fr_IP_para, width = 4, textvariable = textpolar, bg = "white")
	
	labelsnthresh <- tklabel(fr_IP_para, text = " , snthresh: ")
	tk2tip(labelsnthresh, "XCMS: signal to noise ratio cut off")
	textsnthresh <- tclVar(10)
	textsnthreshWidget <- tkentry(fr_IP_para, width = 4, textvariable = textsnthresh, bg = "white")
	
	labelmslevel <- tklabel(fr_IP_para, text = " , mslevel: ")
	tk2tip(labelmslevel, "XCMS: perform peak picking on data given mslevel (n:NULL)")
	textmslevel <- tclVar("n")
	textmslevelWidget <- tkentry(fr_IP_para, width = 6, textvariable = textmslevel, bg = "white")
	
	labelscanrange <- tklabel(fr_IP_para, text = " , scanrange: ")
	tk2tip(labelscanrange, "XCMS: scan range to read (n:NULL)")
	textscanrange <- tclVar("n")
	textscanrangeWidget <- tkentry(fr_IP_para, width = 8, textvariable = textscanrange, bg = "white")
	
	labelnslave <- tklabel(fr_IP_para, text = " , nSlaves: ")
	tk2tip(labelnslave, "number of slaves/cores to be used for parallel peak detection.")
	textnslave <- tclVar(1)
	textnslaveWidget <- ttkcombobox(fr_IP_para, state = "readonly", values = 1:as.numeric(Sys.getenv('NUMBER_OF_PROCESSORS')), width = 4, textvariable = textnslave)
	
	
		
			
			
		
	
	
	tkgrid(labellock, textlockWidget, labelpolar, textpolarWidget, labelsnthresh, textsnthreshWidget, labelmslevel, textmslevelWidget, labelscanrange, textscanrangeWidget,  labelnslave, textnslaveWidget, tklabel(fr_IP_para, text = "    "), sticky = "w")
	
	fr_IP_prof <- tkframe(fr_IP)
	tkgrid(fr_IP_prof, sticky = "w")
	labelprof <- tklabel(fr_IP_prof, text = "                  profmethod: ")
	prof.val = tclVar(1)
	prof.bin <- tkradiobutton(fr_IP_prof, variable = prof.val, value = 1)
	labelprof.bin <- tklabel(fr_IP_prof, text = "bin")
	tk2tip(labelprof.bin, "XCMS: simply bins the intensity into the matrix \ncell closest to it in mass")
	prof.binlin <- tkradiobutton(fr_IP_prof, variable = prof.val, value = 2)
	labelprof.binlin <- tklabel(fr_IP_prof, text = "binlin")
	tk2tip(labelprof.binlin, "XCMS: does the same thing except that it uses\nlinear interpolation to fill in cells that \notherwise would have been left at zero")
	prof.binlinbase <- tkradiobutton(fr_IP_prof, variable = prof.val, value = 3)
	labelprof.binlinbase <- tklabel(fr_IP_prof, text = "binlinbase")
	tk2tip(labelprof.binlinbase, "XCMS: uses linear interpolation between data points")
	prof.intlin <- tkradiobutton(fr_IP_prof, variable = prof.val, value = 4)
	labelprof.intlin <- tklabel(fr_IP_prof, text = "intlin")
	tk2tip(labelprof.intlin, "XCMS: uses integration and linear interpolation \nbetween mass/intensity pairs to determine \nthe equally spaced intensity values")
	tkgrid(labelprof, prof.bin, labelprof.bin, prof.binlin, labelprof.binlin, prof.binlinbase, labelprof.binlinbase, prof.intlin, labelprof.intlin)
	
	
	fr_IP_MF <- tkframe(fr_IP)
	tkgrid(fr_IP_MF, sticky = "w")
	
	IP_method.val <- tclVar(1)
	fr_IP_MF_method <- tkframe(fr_IP_MF)
	tkgrid(fr_IP_MF_method, sticky = "w")
	IP_method.MF <- tkradiobutton(fr_IP_MF_method, variable = IP_method.val, value = 1, command = function(...){
		
		tkconfigure(textfwhmWidget, state = "normal")
		tkconfigure(textMFmaxWidget, state = "normal")
		tkconfigure(textstepWidget, state = "normal")
		tkconfigure(textstepsWidget, state = "normal")
		tkconfigure(textMFmzdiffWidget, state = "normal")
		tkconfigure(labelfwhm, state = "normal")
		tkconfigure(labelMFmax, state = "normal")
		tkconfigure(labelstep, state = "normal")
		tkconfigure(labelsteps, state = "normal")
		tkconfigure(labelMFmzdiff, state = "normal")
		
		tkconfigure(textppmWidget, state = "disable")
		tkconfigure(textpeakwidth.minWidget, state = "disable")
		tkconfigure(textpeakwidth.maxWidget, state = "disable")
		tkconfigure(textprefilter.kWidget, state = "disable")
		tkconfigure(textprefilter.IWidget, state = "disable")
		tkconfigure(textCWmzdiffWidget, state = "disable")
		tkconfigure(textfitgaussWidget, state = "disable")
		tkconfigure(textnoiseWidget, state = "disable")
		tkconfigure(labelppm, state = "disable")
		tkconfigure(labelpeakwidth.min, state = "disable")
		tkconfigure(labelpeakwidth.max, state = "disable")
		tkconfigure(labelprefilter.k, state = "disable")
		tkconfigure(labelprefilter.I, state = "disable")
		tkconfigure(labelCWmzdiff, state = "disable")
		tkconfigure(labelfitgauss, state = "disable")
		tkconfigure(labelnoise, state = "disable")
		
		tkconfigure(PA.MIfilter, state = "normal")
		tkconfigure(labelPA.IAraw, text = "integrated area of original (raw) peak")
		tkconfigure(labelPA.IAfilter, text = "integrated area of filtered peak")
		tkconfigure(labelPA.MIraw, text = "maximum intensity of original (raw) peak")
		tkconfigure(labelPA.MIfilter, text = "maximum intensity of filtered peak", state = "normal")
		
	})
	labelMF <- tklabel(fr_IP_MF_method, text = "MatchedFilter-method")
	tk2tip(labelMF, "XCMS: Find peaks in extracted the chromatographic time domain of the profile matrix.")
	tkgrid(tklabel(fr_IP_MF_method, text = "                "), IP_method.MF, labelMF)
	
	
	fr_IP_MF_p2 <- tkframe(fr_IP_MF)
	tkgrid(fr_IP_MF_p2, sticky = "w")
	
	labelfwhm <- tklabel(fr_IP_MF_p2, text = "                           fwhm: ")
	tk2tip(labelfwhm, "XCMS: full width at half maximum of matched \nfiltration gaussian model peak")
	textfwhm <- tclVar(30)
	textfwhmWidget <- tkentry(fr_IP_MF_p2, width = 4, textvariable = textfwhm, bg = "white")
	
	labelMFmax <- tklabel(fr_IP_MF_p2, text = " , max: ")
	tk2tip(labelMFmax, "XCMS: maximum number of peaks per extracted \nion chromatogram")
	textMFmax <- tclVar(5)
	textMFmaxWidget <- tkentry(fr_IP_MF_p2, width = 4, textvariable = textMFmax, bg = "white")
	
	labelstep <- tklabel(fr_IP_MF_p2, text = " , step: ")
	tk2tip(labelstep, "XCMS: XCMS: step size to use for profile generation")
	textstep <- tclVar(0.1)
	textstepWidget <- tkentry(fr_IP_MF_p2, width = 4, textvariable = textstep, bg = "white")
	
	labelsteps <- tklabel(fr_IP_MF_p2, text = " , steps: ")
	tk2tip(labelsteps, "XCMS: number of steps to merge prior to filtration")
	textsteps <- tclVar(2)
	textstepsWidget <- tkentry(fr_IP_MF_p2, width = 4, textvariable = textsteps, bg = "white")
	
	labelMFmzdiff <- tklabel(fr_IP_MF_p2, text = " , mzdiff: ")
	tk2tip(labelMFmzdiff, "XCMS: minimum difference in m/z for peaks \nwith overlapping retention times")
	textMFmzdiff <- tclVar(0.8)
	textMFmzdiffWidget <- tkentry(fr_IP_MF_p2, width = 4, textvariable = textMFmzdiff, bg = "white")
	
	tkgrid(labelfwhm, textfwhmWidget, labelMFmax, textMFmaxWidget, labelstep, textstepWidget, labelsteps, textstepsWidget, labelMFmzdiff, textMFmzdiffWidget, tklabel(fr_IP_MF_p2, text = "-step*steps          "), sticky = "w")
	

	
	
	
	fr_IP_CW = tkframe(fr_IP)
	tkgrid(fr_IP_CW, sticky = "w")
	
	fr_IP_CW_method = tkframe(fr_IP_CW)
	tkgrid(fr_IP_CW_method, sticky = "w")
	
	IP_method.CW <- tkradiobutton(fr_IP_CW_method, variable = IP_method.val, value = 2, command = function(...){
		tkconfigure(textppmWidget, state = "normal")
		tkconfigure(textpeakwidth.minWidget, state = "normal")
		tkconfigure(textpeakwidth.maxWidget, state = "normal")
		tkconfigure(textprefilter.kWidget, state = "normal")
		tkconfigure(textprefilter.IWidget, state = "normal")
		tkconfigure(textCWmzdiffWidget, state = "normal")
		tkconfigure(textfitgaussWidget, state = "normal")
		tkconfigure(textnoiseWidget, state = "normal")
		tkconfigure(labelppm, state = "normal")
		tkconfigure(labelpeakwidth.min, state = "normal")
		tkconfigure(labelpeakwidth.max, state = "normal")
		tkconfigure(labelprefilter.k, state = "normal")
		tkconfigure(labelprefilter.I, state = "normal")
		tkconfigure(labelCWmzdiff, state = "normal")
		tkconfigure(labelfitgauss, state = "normal")
		tkconfigure(labelnoise, state = "normal")
		
		tkconfigure(textfwhmWidget, state = "disable")
		tkconfigure(textMFmaxWidget, state = "disable")
		tkconfigure(textstepWidget, state = "disable")
		tkconfigure(textstepsWidget, state = "disable")
		tkconfigure(textMFmzdiffWidget, state = "disable")
		tkconfigure(labelfwhm, state = "disable")
		tkconfigure(labelMFmax, state = "disable")
		tkconfigure(labelstep, state = "disable")
		tkconfigure(labelsteps, state = "disable")
		tkconfigure(labelMFmzdiff, state = "disable")
		
		tkconfigure(PA.MIfilter, state = "disable")
		tkconfigure(labelPA.IAraw, text = "integrated peak intensity")
		tkconfigure(labelPA.IAfilter, text = "baseline corrected integrated peak intensity")
		tkconfigure(labelPA.MIraw, text = "maximum peak intensity")
		tkconfigure(labelPA.MIfilter, text = "cannot be chosen", state = "disable")
	})
	labelCW <- tklabel(fr_IP_CW_method, text = "centWave-method")
	tk2tip(labelCW, "XCMS: Peak density and wavelet based feature \ndetection for high resolution LC/MS \ndata in centroid mode.")
	tkgrid(tklabel(fr_IP_CW_method, text = "                "), IP_method.CW, labelCW)
	
	fr_IP_CW_para1 <- tkframe(fr_IP_CW)
	tkgrid(fr_IP_CW_para1, sticky = "w")
	labelppm <- tklabel(fr_IP_CW_para1, text = "                           ppm: ")
	tk2tip(labelppm, "XCMS: maximal tolerated m/z deviation in consecutive \nscans, in ppm (parts per million)")
	textppm <- tclVar("25")
	textppmWidget <- tkentry(fr_IP_CW_para1, width = 4, textvariable = textppm, bg = "white")
	
	labelpeakwidth.min <- tklabel(fr_IP_CW_para1, text = " , peakwidth: min ")
	tk2tip(labelpeakwidth.min, "XCMS: chromatographic peak width, given \nas range (min,max) in seconds")
	labelpeakwidth.max <- tklabel(fr_IP_CW_para1, text = " , max ")
	tk2tip(labelpeakwidth.max, "XCMS: chromatographic peak width, given \nas range (min,max) in seconds")
	textpeakwidth.min <- tclVar("20")
	textpeakwidth.max <- tclVar("50")
	textpeakwidth.minWidget <- tkentry(fr_IP_CW_para1, width = 4, textvariable = textpeakwidth.min, bg = "white")
	textpeakwidth.maxWidget <- tkentry(fr_IP_CW_para1, width = 4, textvariable = textpeakwidth.max, bg = "white")

	
	labelCWmzdiff <- tklabel(fr_IP_CW_para1, text = " , mzdiff: ")
	tk2tip(labelCWmzdiff, "XCMS: minimum difference in m/z for peaks \nwith overlapping retention times, can \nbe negative to allow overlap")
	textCWmzdiff <- tclVar(-0.02)
	textCWmzdiffWidget <- tkentry(fr_IP_CW_para1, width = 5, textvariable = textCWmzdiff, bg = "white")
	
	labelfitgauss <- tklabel(fr_IP_CW_para1, text = " , fitgauss: ")
	tk2tip(labelfitgauss, "XCMS: if T a Gaussian is fitted to each peak")
	textfitgauss <- tclVar("F")
	textfitgaussWidget <- tkentry(fr_IP_CW_para1, width = 4, textvariable = textfitgauss, bg = "white")
	
	labelnoise <- tklabel(fr_IP_CW_para1, text = " , noise: ")
	tk2tip(labelnoise, "XCMS: optional argument which is useful for data that was \ncentroided without any intensity threshold, centroids with \nintensity < noise are omitted from ROI detection")
	textnoise <- tclVar(0)
	textnoiseWidget <- tkentry(fr_IP_CW_para1, width = 4, textvariable = textnoise, bg = "white")
	
	tkgrid(labelppm, textppmWidget, labelpeakwidth.min, textpeakwidth.minWidget, labelpeakwidth.max, textpeakwidth.maxWidget, labelCWmzdiff, textCWmzdiffWidget, labelfitgauss, textfitgaussWidget, labelnoise, textnoiseWidget, tklabel(fr_IP_CW_para1, text = "    "), sticky = "w")
	
	fr_IP_CW_para2 <- tkframe(fr_IP_CW)
	tkgrid(fr_IP_CW_para2, sticky = "w")
	
	labelprefilter.k <- tklabel(fr_IP_CW_para2, text = "                           prefilter: k ")
	tk2tip(labelprefilter.k, "XCMS: Prefilter step for the first phase. Mass traces \nare only retained if they contain at least \n'k' peaks with intensity >= 'I'")
	labelprefilter.I <- tklabel(fr_IP_CW_para2, text = " , I ")
	tk2tip(labelprefilter.I, "XCMS: Prefilter step for the first phase. Mass traces \nare only retained if they contain at least \n'k' peaks with intensity >= 'I'")
	textprefilter.k <- tclVar("3")
	textprefilter.I <- tclVar("100")
	textprefilter.kWidget <- tkentry(fr_IP_CW_para2, width = 4, textvariable = textprefilter.k, bg = "white")
	textprefilter.IWidget <- tkentry(fr_IP_CW_para2, width = 4, textvariable = textprefilter.I, bg = "white")
	tkgrid(labelprefilter.k, textprefilter.kWidget, labelprefilter.I, textprefilter.IWidget, sticky = "w")
	
	tkgrid(tklabel(fr_IP,text="", height = 0, font = fontIntro_para))
	
	tkconfigure(textppmWidget, state = "disable")
	tkconfigure(textpeakwidth.minWidget, state = "disable")
	tkconfigure(textpeakwidth.maxWidget, state = "disable")
	tkconfigure(textprefilter.kWidget, state = "disable")
	tkconfigure(textprefilter.IWidget, state = "disable")
	tkconfigure(textCWmzdiffWidget, state = "disable")
	tkconfigure(textfitgaussWidget, state = "disable")
	tkconfigure(textnoiseWidget, state = "disable")
	
	
	fr_GP <- tkframe(dlg)
	tkpack(fr_GP, side = "top", fill = "x")
	fr_GP_head <- tkframe(fr_GP)
	tkgrid(fr_GP_head, sticky = "w")
	tkgrid(tklabel(fr_GP_head, text = "         Group peaks from different samples (loop)"), sticky = "w")
	
	fr_GP_para <- tkframe(fr_GP)
	tkgrid(fr_GP_para, sticky = "w")
	
	labelbw <- tklabel(fr_GP_para, text = "                  bw: ")
	tk2tip(labelbw, "XCMS: bandwidth (standard deviation or half \nwidth at half maximum) of gaussian \nsmoothing kernel to apply to the peak \ndensity chromatogram")
	textbw <- tclVar("30")
	textbwWidget <- tkentry(fr_GP_para, width = 4, textvariable = textbw, bg = "white")
	
	labelminfrac <- tklabel(fr_GP_para, text = " , minfrac: ")
	tk2tip(labelminfrac, "XCMS: minimum fraction of samples necessary in at \nleast one of the sample groups for it to \nbe a valid group")
	textminfrac <- tclVar("0.5")
	textminfracWidget <- tkentry(fr_GP_para, width = 4, textvariable = textminfrac, bg = "white")
	
	labelminsamp <- tklabel(fr_GP_para, text = " , minsamp: ")
	tk2tip(labelminsamp, "XCMS: minimum number of samples necessary in at \nleast one of the sample groups for it to \nbe a valid group")
	textminsamp <- tclVar("1")
	textminsampWidget <- tkentry(fr_GP_para, width = 4, textvariable = textminsamp, bg = "white")
	
	labelmzwid <- tklabel(fr_GP_para, text = " , mzwid: ")
	tk2tip(labelmzwid, "XCMS: width of overlapping m/z slices to use for \ncreating peak density chromatograms and \ngrouping peaks across samples")
	textmzwid <- tclVar("0.25")
	textmzwidWidget <- tkentry(fr_GP_para, width = 4, textvariable = textmzwid, bg = "white")
	
	labelGPmax <- tklabel(fr_GP_para, text = " , max: ")
	tk2tip(labelGPmax, "XCMS: maximum number of groups to identify in a single m/z slice")
	textGPmax <- tclVar("50")
	textGPmaxWidget <- tkentry(fr_GP_para, width = 4, textvariable = textGPmax, bg = "white")
	
	tkgrid(labelbw, textbwWidget, labelminfrac, textminfracWidget, labelminsamp, textminsampWidget, labelmzwid, textmzwidWidget, labelGPmax, textGPmaxWidget, sticky = "w")
	tkgrid(tklabel(fr_GP_para,text="", height = 0, font = fontIntro_para))
	
	
	
	fr_CR <- tkframe(dlg)
	tkpack(fr_CR, side = "top", fill = "x")
	fr_CR_head <- tkframe(fr_CR)
	tkgrid(fr_CR_head, sticky = "w")
	tkgrid(tklabel(fr_CR_head, text = "         Correct retention time from different samples (loop)"))
	
	fr_CR_para_1 <- tkframe(fr_CR)
	tkgrid(fr_CR_para_1, sticky = "w")
	
	labelmissing <- tklabel(fr_CR_para_1, text = "                  missing: ")
	tk2tip(labelmissing, "XCMS: number of missing samples to allow \nin retention time correction groups")
	textmissing <- tclVar("1")
	textmissingWidget <- tkentry(fr_CR_para_1, width = 4, textvariable = textmissing, bg = "white")
	
	labelextra <- tklabel(fr_CR_para_1, text = " , extra: ")
	tk2tip(labelextra, "XCMS: number of extra peaks to allow in retention \ntime correction correction groups")
	textextra <- tclVar("1")
	textextraWidget <- tkentry(fr_CR_para_1, width = 4, textvariable = textextra, bg = "white")
	
	labelspan <- tklabel(fr_CR_para_1, text = " , span: ")
	tk2tip(labelspan, "XCMS: degree of smoothing for local polynomial regression fitting")
	textspan <- tclVar("0.2")
	textspanWidget <- tkentry(fr_CR_para_1, width = 4, textvariable = textspan, bg = "white")
	
	tkgrid(labelmissing, textmissingWidget, labelextra, textextraWidget, labelspan, textspanWidget, sticky = "w")
	
	fr_CR_para_2 <- tkframe(fr_CR)
	tkgrid(fr_CR_para_2, sticky = "w")
	labelsmooth <- tklabel(fr_CR_para_2, text = "                  smooth: ")
	smooth.val <- tclVar(2)
	smooth.loess <- tkradiobutton(fr_CR_para_2, variable = smooth.val, value = 1)
	smooth.linear <- tkradiobutton(fr_CR_para_2, variable = smooth.val, value = 2)
	labelsmooth.loess <- tklabel(fr_CR_para_2, text = "loess  ")
	tk2tip(labelsmooth.loess, "XCMS: for non-linear alignment")
	labelsmooth.linear <- tklabel(fr_CR_para_2, text = "linear")
	tk2tip(labelsmooth.linear, "XCMS: for linear alignment")
	
	tkgrid(labelsmooth, smooth.loess, labelsmooth.loess, smooth.linear, labelsmooth.linear, sticky = "w")
	
	fr_CR_para_3 <- tkframe(fr_CR)
	tkgrid(fr_CR_para_3, sticky = "w")
	labelfamily <- tklabel(fr_CR_para_3, text = "                  family: ")
	family.val <- tclVar(1)
	family.gaussian <- tkradiobutton(fr_CR_para_3, variable = family.val, value = 1)
	family.symmetric <- tkradiobutton(fr_CR_para_3, variable = family.val, value = 2)
	labelfamily.gaussian <- tklabel(fr_CR_para_3, text = "gaussian  ")
	tk2tip(labelfamily.gaussian, "XCMS: by least-squares with no outlier removal")
	labelfamily.symmetric <- tklabel(fr_CR_para_3, text = "symmetric")
	tk2tip(labelfamily.symmetric, "XCMS: a re-decending M estimator is used with \nTukey's biweight function, allowing \noutlier removal")

	
	tkgrid(labelfamily, family.gaussian, labelfamily.gaussian, family.symmetric, labelfamily.symmetric, sticky = "w")
	tkgrid(tklabel(fr_CR_para_3, text="", height = 0, font = fontIntro_para))
	
	
	
	
	fr_IA <- tkframe(dlg)
	tkpack(fr_IA, fill = "x", side = "top")
	
	fr_IA_head <- tkframe(fr_IA)
	tkgrid(fr_IA_head, sticky = "w")
	tkgrid(tklabel(fr_IA_head, text = "         Integrate areas of missing peaks: "))
	
	fr_IA_para <- tkframe(fr_IA)
	tkgrid(fr_IA_para, sticky = "w")
	
	labelfillpeaks <- tklabel(fr_IA_para, text = "                  fillPeaks: ")
	fillpeaks.val <- tclVar(2)
	fillpeaks.yes <- tkradiobutton(fr_IA_para, variable = fillpeaks.val, value = 1)
	fillpeaks.no <- tkradiobutton(fr_IA_para, variable = fillpeaks.val, value = 2)
	labelfillpeaks.yes <- tklabel(fr_IA_para, text = "Yes  ")
	labelfillpeaks.no <- tklabel(fr_IA_para, text = "No")
	tkgrid(labelfillpeaks, fillpeaks.yes, labelfillpeaks.yes, fillpeaks.no, labelfillpeaks.no, sticky = "w")
	tkgrid(tklabel(fr_IA_para, text="", height = 0, font = fontIntro_para))
	
	
	
	fr_PA <- tkframe(dlg)
	tkpack(fr_PA, side = "top", fill = "x")
	
	fr_PA_head <- tkframe(fr_PA)
	tkgrid(fr_PA_head, sticky = "w")
	tkgrid(tklabel(fr_PA_head, text = "         Peak abundance format: "))
	
	fr_PA_para <- tkframe(fr_PA)
	tkgrid(fr_PA_para, sticky = "w")
	
	PA.val <- tclVar(1)
	PA.IAraw <- tkradiobutton(fr_PA_para, variable = PA.val, value = 1)
	PA.IAfilter <- tkradiobutton(fr_PA_para, variable = PA.val, value = 2)
	PA.MIraw <- tkradiobutton(fr_PA_para, variable = PA.val, value = 3)
	PA.MIfilter <- tkradiobutton(fr_PA_para, variable = PA.val, value = 4)
	labelPA.IAraw <- tklabel(fr_PA_para, text = "integrated area of original (raw) peak")
	labelPA.IAfilter <- tklabel(fr_PA_para, text = "integrated area of filtered peak")
	labelPA.MIraw <- tklabel(fr_PA_para, text = "maximum intensity of original (raw) peak")
	labelPA.MIfilter <- tklabel(fr_PA_para, text = "maximum intensity of filtered peak")
	
	
	
	

	tkgrid(tklabel(fr_PA_para, text = "                  "), PA.IAraw, labelPA.IAraw, PA.MIraw, labelPA.MIraw, sticky = "w")
	tkgrid(tklabel(fr_PA_para, text = "                  "), PA.IAfilter, labelPA.IAfilter, PA.MIfilter, labelPA.MIfilter, sticky = "w")
	tkgrid(tklabel(fr_PA_para, text="", height = 0, font = fontIntro_para), sticky = "w")
	
	
	fr_CAMERA <- tkframe(dlg)
	tkpack(fr_CAMERA, side = "top", fill = "x")

	fr_CAMERA_head <- tkframe(fr_CAMERA)
	tkgrid(fr_CAMERA_head, sticky = "w")
	tkgrid(CAMERA.lab <- tklabel(fr_CAMERA_head, text = "         Peak annotation:"))
	tk2tip(CAMERA.lab, "CAMERA: algorithms for annotation of isotope peaks,\nadducts and fragments in peak lists.")
                        
	fr_CAMERA_para <- tkframe(fr_CAMERA)
	tkgrid(fr_CAMERA_para, sticky = "w")
	CAMERA.iso.val <- tclVar(0)
	CAMERA.frg.val <- tclVar(0)
	CAMERA.add.val <- tclVar(0)
	CAMERA.iso <- tkcheckbutton(fr_CAMERA_para,variable=CAMERA.iso.val)
	CAMERA.frg <- tkcheckbutton(fr_CAMERA_para,variable=CAMERA.frg.val)
	CAMERA.add <- tkcheckbutton(fr_CAMERA_para,variable=CAMERA.add.val)
		CAMERA.add.pn <- tclVar(1)
		CAMERA.add.p <- tkradiobutton(fr_CAMERA_para, variable = CAMERA.add.pn, value = 1)
		CAMERA.add.n <- tkradiobutton(fr_CAMERA_para, variable = CAMERA.add.pn, value = 2)
		labelCAMERA.add.p <- tklabel(fr_CAMERA_para, text="positive")
		labelCAMERA.add.n <- tklabel(fr_CAMERA_para, text="negative")
	labelCAMERA.iso <- tklabel(fr_CAMERA_para, text="isotope")
	labelCAMERA.frg <- tklabel(fr_CAMERA_para, text="fragment")
	labelCAMERA.add <- tklabel(fr_CAMERA_para, text="adduct")

	tkgrid(tklabel(fr_CAMERA_para, text = "                  filter out by selection: "), CAMERA.iso, labelCAMERA.iso, CAMERA.frg, labelCAMERA.frg, CAMERA.add, labelCAMERA.add, tklabel(fr_CAMERA_para, text="( "), CAMERA.add.p, labelCAMERA.add.p, CAMERA.add.n, labelCAMERA.add.n, tklabel(fr_CAMERA_para, text=" )"), sticky = "w")
	tkgrid(tklabel(fr_CAMERA_para, text="", height = 0, font = fontIntro_para), sticky = "w")


	
	onOK <- function()
	{
		ncpu = as.numeric(tclvalue(textnslave))
		if(afterQC){
			if(file.info(tclvalue(textmzdatainput))$isdir){
				
				outpath = paste(tclvalue(textoutput), "/Re-alignment", sep = "")
				dir.create(outpath, showWarnings = F)
				remain_sample = read.table(QCpath, sep = ",", header = T, quote = "\"", as.is = T, skip = 3)
				remain_sample = remain_sample[which(remain_sample[,2]=="No"),1]
				sample_id = gsub("\\d{4,8}_(.*?)|(^[A-MO-Z]\\d{4,8})_(.*?)", "\\1", dir(tclvalue(textmzdatainput), pattern = ".mzXML"))
				sample_id = sub("-", ".", gsub("(.*?).mzXML", "\\1", sample_id))
				if(all(remain_sample %in% sample_id)){
					tkdestroy(dlg)
					tkconfigure(tt,cursor="watch")
					xcmsFiles = dir(tclvalue(textmzdatainput), full.names = T, pattern = ".mzXML")[sort(match(remain_sample, sample_id))]
					pb <<- tkProgressBar("Peak Re-alignment and Annotation - Please wait for alignment processing", "0% done", 0, 100, 0, width = 500)
					cat("Peak Re-alignment and Annotation - Please wait for alignment processing 0 ")
					xset <- switch(tclvalue(IP_method.val), `1` = xcmsSet(xcmsFiles, method = "matchedFilter", profmethod = switch(as.numeric(tclvalue(prof.val)), "bin", "binlin", "binlinbase", "intlin"), lockMassFreq = switch(tclvalue(textlock), "T" = TRUE, "F" = FALSE), mslevel = switch(tclvalue(textmslevel), "n" = NULL, as.numeric(tclvalue(textmslevel))), scanrange = switch(tclvalue(textscanrange), "n" = NULL, as.numeric(tclvalue(textmslevel))), polarity = switch(tclvalue(textpolar), "b" = NULL, "+" = "positive", "-" = "negative"), fwhm = as.numeric(tclvalue(textfwhm)), max = as.numeric(tclvalue(textMFmax)), snthresh = as.numeric(tclvalue(textsnthresh)), step = as.numeric(tclvalue(textstep)), steps = as.numeric(tclvalue(textsteps)), mzdiff = as.numeric(tclvalue(textMFmzdiff))-as.numeric(tclvalue(textstep))*as.numeric(tclvalue(textsteps)), nSlaves = ncpu), `2` = xcmsSet(xcmsFiles, method="centWave", profmethod = switch(as.numeric(tclvalue(prof.val)), "bin", "binlin", "binlinbase", "intlin"), lockMassFreq = switch(tclvalue(textlock), "T" = TRUE, "F" = FALSE), mslevel = switch(tclvalue(textmslevel), "n" = NULL, as.numeric(tclvalue(textmslevel))), scanrange = switch(tclvalue(textscanrange), "n" = NULL, as.numeric(tclvalue(textmslevel))), polarity = switch(tclvalue(textpolar), "b" = NULL, "+" = "positive", "-" = "negative"), ppm = as.numeric(tclvalue(textppm)), peakwidth=as.numeric(c(tclvalue(textpeakwidth.min), tclvalue(textpeakwidth.max))), snthresh = as.numeric(tclvalue(textsnthresh)), prefilter=as.numeric(c(tclvalue(textprefilter.k), tclvalue(textprefilter.I))), mzdiff = as.numeric(tclvalue(textCWmzdiff)), nSlaves = ncpu, fitgauss = switch(tclvalue(textfitgauss), "T" = TRUE, "F" = FALSE), noise = as.numeric(tclvalue(textnoise))))
					
					setTkProgressBar(pb, value = 100, "Peak Re-alignment and Annotation - Please wait for alignment processing (100% done)", "Finished 100%")
					Sys.sleep(1)
					close(pb)
					textoutput <<- tclVar(tclvalue(textoutput))
					align_loop(loop = 1, xset = xset, outpath = outpath, bw = as.numeric(tclvalue(textbw)), minfrac = as.numeric(tclvalue(textminfrac)), minsamp = as.numeric(tclvalue(textminsamp)), mzwid = as.numeric(tclvalue(textmzwid)), max = as.numeric(tclvalue(textGPmax)), missing_n = as.numeric(tclvalue(textmissing)), extra = as.numeric(tclvalue(textextra)), smooth = ifelse(tclvalue(smooth.val)==1, "loess", "linear"), span = as.numeric(tclvalue(textspan)), family = ifelse(tclvalue(family.val)==1, "gaussian", "symmetric"), fillpeak = (tclvalue(fillpeaks.val)==1), PAformat = as.numeric(tclvalue(PA.val)), IP_method = as.numeric(tclvalue(IP_method.val)), filter.iso = as.numeric(tclvalue(CAMERA.iso.val)), filter.frg = as.numeric(tclvalue(CAMERA.frg.val)), filter.add = as.numeric(tclvalue(CAMERA.add.val)), add.pn = as.numeric(tclvalue(CAMERA.add.pn)))
					text_missing <<- "NA"
				}else{
					tkmessageBox(title = "Error", message = "Some IDs' mzfile is not found. Please check all IDs' mzfile are in the input folder.", icon = "error", type = "ok")
					cat("Peak Re-alignment and Annotation(Error) - Some IDs' mzfile is not found. Please check all IDs' mzfile are in the input folder.\n")
					tkfocus(dlg)
				}
			}else{
				tkmessageBox(title = "Error", message = "Please check the type of all files in the input folder is in mzXML format.", icon = "error", type = "ok")
				cat("Peak Re-alignment and Annotation(Error) - Please check the type of all files in the input folder is in mzXML format.\n")
				tkfocus(dlg)
			}
		}else{
			outpath = paste(tclvalue(textoutput), "/Peak Alignment and Annotation", sep = "")
			dir.create(outpath, showWarnings = F)
			tkdestroy(dlg)
			xcmsFiles = dir(tclvalue(textmzdatainput), full.names = T, pattern = ".mzXML")
			pb <<- tkProgressBar("Peak Alignment and Annotation - Please wait for alignment processing", "0% done", 0, 100, 0, width = 500)
			cat("Peak Alignment and Annotation - Please wait for alignment processing \n0 ")
			xset <- switch(tclvalue(IP_method.val), `1` = xcmsSet(xcmsFiles, method = "matchedFilter", profmethod = switch(as.numeric(tclvalue(prof.val)), "bin", "binlin", "binlinbase", "intlin"), lockMassFreq = switch(tclvalue(textlock), "T" = TRUE, "F" = FALSE), mslevel = switch(tclvalue(textmslevel), "n" = NULL, as.numeric(tclvalue(textmslevel))), scanrange = switch(tclvalue(textscanrange), "n" = NULL, as.numeric(tclvalue(textmslevel))), polarity = switch(tclvalue(textpolar), "b" = NULL, "+" = "positive", "-" = "negative"), fwhm = as.numeric(tclvalue(textfwhm)), max = as.numeric(tclvalue(textMFmax)), snthresh = as.numeric(tclvalue(textsnthresh)), step = as.numeric(tclvalue(textstep)), steps = as.numeric(tclvalue(textsteps)), mzdiff = as.numeric(tclvalue(textMFmzdiff))-as.numeric(tclvalue(textstep))*as.numeric(tclvalue(textsteps)), nSlaves = ncpu), `2` = xcmsSet(xcmsFiles, method="centWave", profmethod = switch(as.numeric(tclvalue(prof.val)), "bin", "binlin", "binlinbase", "intlin"), lockMassFreq = switch(tclvalue(textlock), "T" = TRUE, "F" = FALSE), mslevel = switch(tclvalue(textmslevel), "n" = NULL, as.numeric(tclvalue(textmslevel))), scanrange = switch(tclvalue(textscanrange), "n" = NULL, as.numeric(tclvalue(textmslevel))), polarity = switch(tclvalue(textpolar), "b" = NULL, "+" = "positive", "-" = "negative"), ppm = as.numeric(tclvalue(textppm)), peakwidth = as.numeric(c(tclvalue(textpeakwidth.min), tclvalue(textpeakwidth.max))), snthresh=as.numeric(tclvalue(textsnthresh)), prefilter=as.numeric(c(tclvalue(textprefilter.k), tclvalue(textprefilter.I))), mzdiff = as.numeric(tclvalue(textCWmzdiff)), nSlaves = ncpu, fitgauss = switch(tclvalue(textfitgauss), "T" = TRUE, "F" = FALSE), noise = as.numeric(tclvalue(textnoise))))
			
			setTkProgressBar(pb, value = 100, "Peak Alignment and Annotation - Please wait for alignment processing (100% done)", "Finished 100%")
			Sys.sleep(1)
			close(pb)
			align_loop(loop = 1, xset = xset, outpath = outpath, bw = as.numeric(tclvalue(textbw)), minfrac = as.numeric(tclvalue(textminfrac)), minsamp = as.numeric(tclvalue(textminsamp)), mzwid = as.numeric(tclvalue(textmzwid)), max = as.numeric(tclvalue(textGPmax)), missing_n = as.numeric(tclvalue(textmissing)), extra = as.numeric(tclvalue(textextra)), smooth = ifelse(tclvalue(smooth.val)==1, "loess", "linear"), span = as.numeric(tclvalue(textspan)), family = ifelse(tclvalue(family.val)==1, "gaussian", "symmetric"), fillpeak = (tclvalue(fillpeaks.val)==1), PAformat = as.numeric(tclvalue(PA.val)), IP_method = as.numeric(tclvalue(IP_method.val)), filter.iso = as.numeric(tclvalue(CAMERA.iso.val)), filter.frg = as.numeric(tclvalue(CAMERA.frg.val)), filter.add = as.numeric(tclvalue(CAMERA.add.val)), add.pn = as.numeric(tclvalue(CAMERA.add.pn)))
		}
		tkconfigure(tt,cursor="arrow")
		tkmessageBox(title=ALtitle,message="RT and M/Z correcting.", icon = "info", type = "ok")
		cat(ALtitle,  " - RT and M/Z correcting.\n", sep = "")
	}
	onCancel <- function()
	{
		
		
		tkdestroy(dlg)
		tkfocus(tt)
	}
	fr <- tkframe(dlg)
	tkpack(fr, side = "top")
	OK.but     <-tkbutton(fr,text="     OK     ",command=onOK)
	Cancel.but <-tkbutton(fr,text="   Cancel   ",command=onCancel)
	tkgrid(tklabel(fr,text="               "), OK.but,tklabel(fr,text="                                                         "), Cancel.but, tklabel(fr,text="               "))
	tkgrid(tklabel(fr,text="    ", height = 0, font = fontIntro_para))

	tkfocus(xcmscitation)
	tkbind(dlg, "<Destroy>", function() {tkgrab.release(dlg);tkfocus(tt)})
	
	
	tkwait.window(dlg)

	
}  


pidBatch <- function(){

	outpath <- paste(ifelse(exists("textoutput"),tclvalue(textoutput),getwd()), "/Peak Identification", sep = "")
	dir.create(outpath, showWarnings = F)
		
    load("dbMetabo.RData")
	

	pidSearch <- function(){
	 string <- gsub("\n{2,}","\n",gsub("[a-df-zA-DF-Z!\"#$%&'()*+,/:;<=>?@\\\\^_`{|}~]|\\[|\\]|0.01 ","",tclvalue(tkget(mass.txt,"0.0","end"))))
	 tkconfigure(mass.txt,state="normal",endline=1)
	 tkinsert(mass.txt,"end",string)
	 tkconfigure(mass.txt,state="normal",endline=str_count(string,"\n")+1)

	 if(tclvalue(mass.in)!=""){
	  sep <- ifelse(grepl("\t",readLines(tclvalue(mass.in),n=1)),"\t",",")
	  peaks <- read.table(tclvalue(mass.in),sep=sep,header=T,stringsAsFactors=F)
	  target <- peaks[,ifelse(ncol(peaks)>1,2,1)]
	 } else {
	  target <- sapply(strsplit(tclvalue(tkget(mass.txt,"0.0","end")),"\n"),as.numeric)
	   target <- target[!is.na(target)]
	 }

	
	
	

	 adduct <- switch(tclvalue(ion.choose),N="","+"=ion.pos,"-"=ion.neg)[as.numeric(tkcurselection(ion.lst))+1]
	
	
		
	 result <- sapply(target,function(tM)
	 {
	  tol <- as.numeric(tclvalue(tol.val))
	
	  tol <- ifelse(as.character(tclvalue(tol.type))=="Da",tol,tM*tol/1e6)

	  if(length(adduct)>0){
	   do.call("rbind",sapply(adduct,function(a){

		 M <- HMDB$Mass; Mass <- eval(parse(text=sub("([0-9]+)M","\\1 \\* M",Adduct$IonMass[Adduct$Mode==a])))
		 abs.delta <- abs(Mass-tM)
		 idx <- which(abs.delta<=tol)
		 hmdb <- HMDB[idx,]
		  if(nrow(hmdb)>0) hmdb <- data.frame(DB="HMDB",hmdb,abs.delta=round(abs.delta[idx],10),stringsAsFactors=F)

		 M <- MB$Mass; Mass <- eval(parse(text=sub("([0-9]+)M","\\1 \\* M",Adduct$IonMass[Adduct$Mode==a])))
		 abs.delta <- abs(Mass-tM)
		 idx <- which(abs.delta<=tol & MB$Mode%in%ifelse(tclvalue(ion.choose)=="N","",tclvalue(ion.choose)))
		 mb <- MB[idx,c("ID","Name","Chem.Formu","Mass")]
		  if(nrow(mb)>0) mb <- data.frame(DB="MassBank",mb,abs.delta=round(abs.delta[idx],10),stringsAsFactors=F)

		 tmp <- rbind(hmdb,mb)

		 if(nrow(tmp)>0) data.frame(Target=tM,Adduct=a,tmp,stringsAsFactors=F)
		 else NULL
		},simplify=F)
	   )
	  } else {
	   abs.delta <- abs(HMDB$Mass-tM)
	   hmdb <- HMDB[which(abs.delta<=tol),]
	   if(nrow(hmdb)>0) data.frame(Target=tM,Adduct.Model="M",DB="HMDB",hmdb,abs.delta=abs.delta[which(abs.delta<=tol)],stringsAsFactors=F)
	   else NULL
	  }
	 },simplify=F)
	 result <- do.call("rbind",result)


	 if(length(result)>0){
	  if(tclvalue(mass.in)!=""){
	   write.table(cbind(peaks[match(result[,1],peaks[,2]),1:3],result[,-1]),sprintf("%s/peakID_search_results_%s.csv",outpath,gsub(":","-",as.character(Sys.time()))),sep=",",row.name=F)
	  } else {
	   write.table(result,sprintf("%s/peakID_search_results_%s.csv",outpath,gsub(":","-",as.character(Sys.time()))),sep=",",row.name=F)
	  }

	  result$abs.delta <- format(result$abs.delta,digit=5,scientific=T)

	  db <- tktoplevel(); if(isIcon) tk2ico.set(db,icon)
	  tktitle(db) <- "Peak Identification"

	  list_array <- tclArray()
	  for(j in 1:ncol(result))
	   list_array[[0,j-1]] <- as.tclObj(c("Target","Adduct","Database","ID","Name","Chemical Formula","Mass","| Delta |")[j],drop=T)
	  for(i in 1:nrow(result))
	   for(j in 1:ncol(result))
		list_array[[i,j-1]] <- as.tclObj(result[i,j],drop=T)

	  tab.vscr <- tkscrollbar(db,repeatinterval=20,orient="vertical",command=function(...) tkyview(table_list,...))
	  tab.hscr <- tkscrollbar(db,repeatinterval=1,orient="horizontal",command = function(...) tkxview(table_list,...))

      table_list <- tk2table(db,variable=list_array,rows=nrow(result)+1,cols=ncol(result),titlerows=1,background="white",yscrollcommand=function(...) tkset(tab.vscr, ...),xscrollcommand=function(...) tkset(tab.hscr, ...),resizeborders="col",bordercursor="sb_h_double_arrow",padx=2,height=15,maxwidth=1e6)

	  tcl(table_list,"tag","configure","coltitle",anchor='w')
	  tkgrid(table_list,tab.vscr)
	  tkgrid.configure(tab.vscr,sticky="nsw")
	  tkgrid(tab.hscr,sticky="new")

	  nChar <- apply(result,2,function(x) max(nchar(x)))
	  sapply(1:length(nChar),function(x) tcl(.Tk.ID(table_list),"width",as.character(x-1),as.character(nChar[x]+10)))
	

	  colWidth <- tclvalue(tcl(.Tk.ID(table_list),"width"))

	  tkconfigure(table_list,selectmode="extended",rowseparator="\"\n\"",colseparator="\"\t\"")
	  tkbind(table_list,"<ButtonRelease-1>",function(){
	   if(tclvalue(tcl(.Tk.ID(table_list),"width"))==colWidth){
		sel <- sapply(strsplit(tclvalue(tkcurselection(table_list)),","),as.integer)
		if(length(sel)==2 & sel[2]==3){
		 switch(result[sel[1],"DB"],HMDB=shell(sprintf("start http://www.hmdb.ca/metabolites/%s",result[sel[1],"ID"])),
									MassBank=shell(sprintf("start http://www.massbank.jp/jsp/Dispatcher.jsp?type=disp\"&\"id=%s\"&\"site=%s",result[sel[1],"ID"],MB$Site[MB$ID==result[sel[1],"ID"]])))
		}
	   }
	   colWidth <<- tclvalue(tcl(.Tk.ID(table_list),"width"))
	  })
	 } else {
	  tkmessageBox(title="Information",message="No matched results",icon="info",type="ok")
	 }
	}

	pid.invoke <- F

	dlg <- tktoplevel(); if(isIcon) tk2ico.set(dlg,icon)
	tkwm.title(dlg, "Batch Search")
	fr_basic <- tkframe(dlg)

	check <- function(){
		all(	gsub("\n","",tclvalue(tkget(mass.txt,"0.0","end")))!="",
			ifelse(tclvalue(ion.choose)=="N",T,tclvalue(tkcurselection(ion.lst))!=""),
			suppressWarnings(!is.na(as.numeric(tclvalue(tkget(tol.box)))))	)
	}

	tkgrid(tklabel(fr_basic,text="",height=1))
	fr_mass <- tkframe(fr_basic)
		fr_mass_1 <- tkframe(fr_mass)
			mass.in <- tclVar("")
			mass.lab <- tklabel(fr_mass_1,text="Mass (Da):",width=10)
			mass.box <- tkentry(fr_mass_1,width=40,bg="white",textvariable=mass.in,xscrollcommand=function(...) {
				tkconfigure(mass.txt,state=ifelse(tclvalue(mass.in)=="","normal","disable"))
				
			})
			mass.but <- tkbutton(fr_mass_1,text="...",command=function() tclvalue(mass.in) <- tkgetOpenFile(initialdir=tclvalue(mass.in),filetypes = "{{Text Files} {.txt .csv}}"))
				tkgrid(mass.lab,mass.box,mass.but)
				tkgrid.configure(mass.lab,mass.box,sticky="w")
		fr_mass_2 <- tkframe(fr_mass)
			mass.vscr <- tkscrollbar(fr_mass_2,repeatinterval=5,command=function(...)tkyview(mass.txt,...))
			mass.hscr <- tkscrollbar(fr_mass_2,repeatinterval=5,orient="horizontal",command=function(...)tkxview(mass.txt,...))
			mass.txt <- tktext(fr_mass_2,width=30,height=10,state="normal",xscrollcommand=function(...)tkset(mass.hscr,...),yscrollcommand=function(...)tkset(mass.vscr,...),bg="white")
			tkinsert(mass.txt,"end","176.118 ") 
			
			
				
				
				
			

				tkgrid(mass.txt,mass.vscr)
				tkgrid(mass.hscr)
				tkgrid.configure(mass.txt,sticky="w")
				tkgrid.configure(mass.vscr,sticky="ns")
				tkgrid.configure(mass.hscr,sticky="ew")
		tkgrid(fr_mass_1)
		tkgrid(fr_mass_2)
		tkgrid.configure(fr_mass_1,fr_mass_2,sticky="ns")
	tkgrid(fr_mass)
	tkgrid.configure(fr_mass,sticky="w")

	tkgrid(tklabel(fr_basic,text="",height=1))
	fr_ion <- tkframe(fr_basic)
		ion.pos <- sort(c("M+3H","M+2H+Na","M+H+2Na","M+3Na","M+2H","M+H+NH4","M+H+Na","M+H+K","M+ACN+2H","M+2Na","M+2ACN+2H","M+3ACN+2H","M+H","M+NH4","M+Na","M+CH3OH+H","M+K","M+ACN+H","M+2Na+H","M+IsoProp+H","M+ACN+Na","M+2K+H","M+DMSO+H","M+2ACN+H","M+2Na-H","M+IsoProp+Na+H","2M+H","2M+NH4","2M+Na","2M+3H2O+2H","2M+K","2M+ACN+H","2M+ACN+Na"),decreasing=T)
		ion.neg <- sort(c("M-3H","M-2H","M-H2O-H","M-H","M+Na-2H","M+Cl","M+K-2H","M+FA-H","M+Hac-H","M+Br","M+TFA-H","2M-H","2M+FA-H","2M+Hac-H","3M-H"),decreasing=T)
		ion.choose <- tclVar("+") 
		fr_ion_1 <- tkframe(fr_ion)
			ion <- c("Neutral","Positive","Negative")
			ion.lab <- tklabel(fr_ion_1,text="Ion Mode:",width=10)
			ion.rb1.lab <- tklabel(fr_ion_1,text="Neutral")
			ion.rb2.lab <- tklabel(fr_ion_1,text="Positive")
			ion.rb3.lab <- tklabel(fr_ion_1,text="Negative")
			ion.rb1 <- tkradiobutton(fr_ion_1,command=function() {tclvalue(ion.ls) <- ""; })
			ion.rb2 <- tkradiobutton(fr_ion_1,command=function() {tclvalue(ion.ls) <- ion.pos; tkselection.clear(ion.lst,"0","100"); tkselection.set(ion.lst,9)})
			ion.rb3 <- tkradiobutton(fr_ion_1,command=function() {tclvalue(ion.ls) <- ion.neg; tkselection.clear(ion.lst,"0","100"); tkselection.set(ion.lst,8)})
				tkconfigure(ion.rb1,variable=ion.choose,value="N")
				tkconfigure(ion.rb2,variable=ion.choose,value="+")
				tkconfigure(ion.rb3,variable=ion.choose,value="-")
				tkgrid(ion.lab,ion.rb1,ion.rb1.lab,ion.rb2,ion.rb2.lab,ion.rb3,ion.rb3.lab)
				tkgrid.configure(ion.lab,ion.rb1,ion.rb1.lab,ion.rb2,ion.rb2.lab,ion.rb3,ion.rb3.lab,sticky="w")
		fr_ion_2 <- tkframe(fr_ion)
		ion.ls <- tclVar(ion.pos) 
		ion.lab2 <- tklabel(fr_ion_2,text="Adduct:    ",width=10)
		ion.vscr <- tkscrollbar(fr_ion_2,repeatinterval=5,command=function(...)tkyview(ion.lst,...))
		ion.hscr <- tkscrollbar(fr_ion_2,repeatinterval=5,orient="horizontal",command=function(...)tkxview(ion.lst,...))
		ion.lst <- tklistbox(fr_ion_2,width=18,height=10,selectmode="extended",listvariable=ion.ls,xscrollcommand=function(...) tkset(ion.hscr,...),yscrollcommand=function(...)tkset(ion.vscr,...),background="white",cursor="hand2")
		tkselection.set(ion.lst,9) 
			selectAll <- function() tkselection.set(ion.lst,"0","100")
			tkbind(ion.lst,"<Control-A>",selectAll)
			tkbind(ion.lst,"<Control-a>",selectAll)
			tkbind(ion.lst,"<ButtonRelease-1>",function() { tkconfigure(OK.but,state=ifelse(check(),"normal","disable")) })
			tkgrid(ion.lab2,ion.lst,ion.vscr)
			tkgrid(tklabel(fr_ion_2,text="",width=10),ion.hscr)
			tkgrid.configure(ion.lab2,sticky="n")
			tkgrid.configure(ion.vscr,sticky="ns")
			tkgrid.configure(ion.hscr,sticky="ew")
		tkgrid(fr_ion_1)
		tkgrid(fr_ion_2)
		tkgrid.configure(fr_ion_1,fr_ion_2,sticky="w")
	tkgrid(fr_ion)
	tkgrid.configure(fr_ion,sticky="w")

	tkgrid(tklabel(fr_basic,text="",height=1))
	fr_tol <- tkframe(fr_basic)
		tol.val <- tclVar("0.01 ")

		tol.type <- tclVar("Da")
		tol.lab <- tklabel(fr_tol,text="Tolerance:",width=10)
		tol.box <- tkentry(fr_tol,width=20,textvariable=tol.val,bg="white",xscrollcommand=function(...){
		 
			string <- tclvalue(tkget(tol.box))
			string <- ifelse(string=="0.01 ",string,gsub("[a-df-zA-DF-Z!\"#$%&'()*+,/:;<=>?@\\\\^_`{|}~]|\\[|\\]","",string))
			tclvalue(tol.val) <- string

			
		 
		})
		
			
			
		

		tol.rb1.lab <- tklabel(fr_tol,text="Da")
		tol.rb2.lab <- tklabel(fr_tol,text="ppm")
		tol.rb1 <- tkradiobutton(fr_tol)
		tol.rb2 <- tkradiobutton(fr_tol)
			tkconfigure(tol.rb1,variable=tol.type,value="Da")
			tkconfigure(tol.rb2,variable=tol.type,value="ppm")

			tkgrid(tol.lab,tol.box,tklabel(fr_tol,text="",width=1),tol.rb1,tol.rb1.lab,tol.rb2,tol.rb2.lab)

			tkgrid.configure(tol.lab,tol.box,sticky="ns")
	tkgrid(fr_tol)
	tkgrid.configure(fr_tol,sticky="w")
	tkgrid(tklabel(fr_basic,text="",height=1))

	fr_detail <- tkframe(dlg)



	
	tkgrid(fr_basic,fr_detail)

	tkgrid(tklabel(dlg,text="",height=1))
	fr <- tkframe(dlg)
	Reset.but <- tkbutton(fr,text="  Reset   ",command=function(){
		tclvalue(ion.choose) <- "N"
		tclvalue(ion.ls) <- ""
		tclvalue(tol.val) <- ""
		tclvalue(tol.type) <- "Da"
	})
	OK.but <- tkbutton(fr,text="  Search  ",command=pidSearch)
	
	Cancel.but <- tkbutton(fr,text="  Cancel  ",command=function(){ tkdestroy(dlg); tkfocus(tt)})
	tkgrid(tklabel(fr,text="   "),Reset.but,tklabel(fr,text="                      "),OK.but,tklabel(fr,text="                      "), Cancel.but, tklabel(fr,text="   "))
	tkgrid(tklabel(fr,text="    ", font=fontIntro_para, height=0))
	tkgrid(fr)

	tkgrid(tklabel(dlg,text="",height=1))
	tkfocus(dlg)
	tkbind(dlg, "<Destroy>", function() {tkgrab.release(dlg);tkfocus(tt)})

	pid.invoke <- F
	tkwait.window(dlg)
}


