
rm(list = ls())
gc()
library(tcltk)

cat("OS:",Sys.info()["sysname"],Sys.info()["release"],"\n")
cat("R:",R.Version()$version.string,R.Version()$arch,"\n")
cat("Tcl:\t",Sys.getenv("tclDir"),"\n")

 
 






libname = .packages(all.available = TRUE)




###
# R<3.5.0: biocLite, else BiocManager
ip <- ifelse(.Platform$OS.type=="windows",system("ping -n 1 cran.r-project.org"),system("ping -c 1 cran.r-project.org"))
if(ip==0){
  if(as.numeric(R.Version()$major)<=3&&as.numeric(R.Version()$minor)<5){
    source("http://bioconductor.org/biocLite.R")
  }
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
}
###


osIsWin <- Sys.info()["sysname"]=="Windows"
if(osIsWin & as.numeric(R.Version()$major)>=3 & as.numeric(R.Version()$minor)>=4.0){ # R-3.4.0 or later
 library(tcltk)
 addTclPath(paste(gsub("\\\\","/",Sys.getenv("tclDir")),"lib",sep="/")) # addTclPath("C:/ActiveTcl/lib")
} else {
 Sys.setenv(MY_TCLTK="yes")
 library(tcltk)
}

if(any(!c('tcltk', 'tcltk2', 'limma', 'mzR', 'xcms', 'tkrplot', 'graphics', 'scatterplot3d', 'car', 'png', 'fields', 'boot', 'snow', 'gtools','data.table') %in% libname) & ip==1){
	tkmessageBox(title = "Error", message = "Not able to connect to Internet.\nPlease check your network connection device.", icon = "error", type = "ok")
	cat("MainInterface(Error) - Not able to connect to Internet.\nPlease check your network connection device.\n")
	stop()
}

###
pkg <- c("limma","mzR","xcms","CAMERA","tcltk2","tcltk","tkrplot","graphics","scatterplot3d","car","png","fields","boot",
"grid","parallel","doSNOW","foreach","gtools","ropls","matrixStats","stringr","dplyr","SPIA",'data.table',"Hmisc","sjmisc","tkImgR","investr","tm","wordcloud")

# if(!"devtools" %in% libname) install.packages("devtools", repos="http://cran.csie.ntu.edu.tw")
# if(!"remotes" %in% libname) install.packages("remotes", repos="http://cran.csie.ntu.edu.tw")
# # library(devtools)
# library(remotes)
# install_version("rlang",version="0.4.10")


pkg_bio <- c("limma","mzR","xcms","CAMERA","ropls","matrixStats","stringr","SPIA")

pkg_cran <- setdiff(pkg,c(pkg_bio,"tcltk","grid"))
pkg_cran <- pkg_cran[!pkg_cran%in%libname]

if(length(pkg_cran)!=0){install.packages(pkg_cran, repos = "http://cran.csie.ntu.edu.tw", dependencies = TRUE)}

pkg_bio <- pkg_bio[!pkg_bio%in%libname]
if(length(pkg_bio)!=0){if(as.numeric(R.Version()$major)<=3 && as.numeric(R.Version()$minor)<5){
  biocLite(pkg_bio,dep=T, ask = F)
}else{
  BiocManager::install(pkg_bio, dep=T, ask = F)
}
}

sapply(pkg, require, character.only = TRUE)
tclRequire("BWidget")
tclRequire("Tktable")
tclRequire("Img")
###


library(limma)  
library(mzR)
library(xcms)
library(CAMERA)
library(tcltk2)
library(tcltk)
tclRequire("BWidget")
tclRequire("Tktable")
tclRequire("Img")
library(tkrplot)
library(graphics)
library(scatterplot3d)
library(car)
library(png)
library(fields)
library(boot)
library(grid)
library(parallel)
library(doSNOW)
library(foreach)

library(gtools)
library(ropls)
library(matrixStats)
library(stringr)

groups <- xcms:::groups	

mass <- NULL 

	fontHeading = tkfont.create(family = "times", size = 14, weight = "bold", slant = "italic")
	fontIntro = tkfont.create(family = "arial", size = 10)
	fontIntro_para = tkfont.create(family = "arial", size = 1)
	fontsubtitle = tkfont.create(family = "times", size = 12, weight = "bold")
	fontTextLabel = tkfont.create(family = "times", size = 10)
	fontTextLabel_para = tkfont.create(family = "times", size = 11)
	fontRun = tkfont.create(family = "courier", size = 10, weight = "bold", underline = F, overstrike = F)

# print(Sys.getenv("subfunc"))
#source(paste("W:/SMART/Program", "Nagisa.r", sep = "/"))
source("SMART_V1.2_Sub.r")







tt <- tktoplevel()
tkwm.title(tt, "SMART")

if(osIsWin){
	icon <- tk2ico.create(paste("Ushio.ico", sep = "/"))
	tk2ico.set(tt,icon)
}
	
	
ws = which(c(as.numeric(tkwm.maxsize(tt))[2]>=900, (as.numeric(tkwm.maxsize(tt))[2]<900 & as.numeric(tkwm.maxsize(tt))[2]>=768), as.numeric(tkwm.maxsize(tt))[2]<768))


tkpack(listframe <- tkframe(tt,width = switch(ws, 210, 170, 170), height = switch(ws, 600, 480, 420)), side = "left")

tkpack(plotframe <- tkframe(tt,width = switch(ws, 800, 760, 720), height = switch(ws, 600, 480, 420)), side = "left")
yScr <- tkscrollbar(listframe, command = function(...)tkyview(treeWidget,...))
xScr <- tkscrollbar(listframe, command = function(...)tkxview(treeWidget,...), orient="horizontal")

treeWidget <- tkwidget(listframe, "Tree", yscrollcommand = function(...)tkset(yScr,...), xscrollcommand = function(...)tkset(xScr,...), width = switch(ws, 25, 20, 20), height = switch(ws, 51, 43, 40))

tkgrid(treeWidget,yScr)
tkgrid.configure(yScr, sticky="nsw")
tkgrid(xScr, sticky = "new")






topMenu <- tkmenu(tt)              
tkconfigure(tt, menu = topMenu)    
fileMenu <- tkmenu(topMenu, tearoff = FALSE)
viewMenu <- tkmenu(topMenu, tearoff = FALSE)

peakMenu <- tkmenu(topMenu, tearoff = FALSE)
	AlignsubMenu <- tkmenu(topMenu, tearoff = FALSE)  
	DatasubMenu <- tkmenu(topMenu, tearoff = FALSE)  
	QCsubMenu <- tkmenu(topMenu, tearoff = FALSE)  
statMenu <- tkmenu(topMenu, tearoff = FALSE)
	batchsubMenu <- tkmenu(topMenu, tearoff = FALSE)  
peakInfoMenu <- tkmenu(topMenu, tearoff=F)


tkadd(topMenu, "cascade", label = "File", menu = fileMenu) 
tkadd(fileMenu, "command", label = "Data Import", command = function(){
	DataImport()
})



tkadd(topMenu, "cascade", label = "Viewer", menu = viewMenu)

tkadd(viewMenu, "command", label="2D Plot", command = function(){
	if(length(unlist(strsplit(tclvalue(tcl(treeWidget,"selection","get")), split = " |,,")))==1){
		viewer <<- "2D"
		tkconfigure(save.but, state = "normal")
		eval(show2D())
	}
})
tkadd(viewMenu, "command", label="3D Plot", command = function(){
	if(length(unlist(strsplit(tclvalue(tcl(treeWidget,"selection","get")), split = " |,,")))==1){
		viewer <<- "3D"
		tkconfigure(save.but, state = "normal")
		eval(show3D())
	}
})
tkadd(viewMenu, "command", label="TIC Plot", command = function(){
	if(length(unlist(strsplit(tclvalue(tcl(treeWidget,"selection","get")), split = " |,,")))){
		viewer <<- "TIC"
		tkconfigure(save.but, state = "normal")
		eval(showTIC())
	}
})

tkadd(viewMenu, "command", label = "TIC Clustering", command = function(){
	TICclustering()
})

tkentryconfigure(viewMenu, 0, state = "disable")
tkentryconfigure(viewMenu, 1, state = "disable")
tkentryconfigure(viewMenu, 2, state = "disable")
tkentryconfigure(viewMenu, 3, state = "disable")


tkadd(topMenu, "cascade", label = "Peak Analysis", menu = peakMenu)


tkadd(peakMenu, "command", label = "Peak Alignment and Annotation", command = function(){
	Align(ALtitle = "Peak Alignment and Annotation")
})



tkadd(peakMenu, "command", label = "Data Preprocessing", command = function(){
	DataPrepro()
})



tkadd(peakMenu, "cascade", label = "Quality Control", menu = QCsubMenu)
	
	tkadd(QCsubMenu, "command", label = "Peak/Sample filtering", command = function(){
		QualityControl()
	})

	
	tkadd(QCsubMenu, "command", label = "Output Viewer", command = function(){
		QualityCtrlViewer()
	})

tkentryconfigure(peakMenu, 0, state = "disable")
tkentryconfigure(peakMenu, 1, state = "disable")
tkentryconfigure(QCsubMenu, 0, state = "disable")




tkadd(topMenu, "cascade", label = "Statistical Methods", menu = statMenu)

tkadd(statMenu, "cascade", label = "Batch Effect Detection", menu = batchsubMenu)

tkadd(batchsubMenu, "command", label = "Principal Component Analysis (PCA)", command = function(){
	pcaStat()
})



tkadd(batchsubMenu, "command", label = "Latent Group (LG)", command = function(){
	LatentGroup()
})



tkadd(statMenu, "command", label = "Analysis of Covariance (ANCOVA)", command = function(){
	ANCOVA()
})



tkadd(statMenu, "command", label = "Partial Least Squares (PLS/PLS-DA)", command = function(){
	PLS()
})

tkentryconfigure(batchsubMenu, 0, state = "disable")
tkentryconfigure(batchsubMenu, 1, state = "disable")
tkentryconfigure(statMenu, 1, state = "disable")
tkentryconfigure(statMenu, 2, state = "disable")



tkadd(topMenu, "cascade", label = "Peak Identification", menu = peakInfoMenu) 
tkadd(peakInfoMenu, "command", label = "Batch Search", command = function(){
	pidBatch()
})



m_reso <- 100
t_reso <- 100
lb_inten <- -Inf
ub_inten <- Inf
rt_begin <- 20
rt_range <- 50
mz_begin <- 100
mz_range <- 500
h_angle <- 15
v_angle <- 15

tkpack(tfr <- tkframe(plotframe), side = "bottom", fill = "x")
tkpack(tklabel(tfr, text = "", height = 0, font = fontIntro_para), side = "bottom")

refresh.but <- tkbutton(tfr, text = "Refresh", width = 18)
tkpack(refresh.but, side = "left", anchor = "s")
tkpack(tklabel(tfr, text = "      ",width = switch(ws, 27, 25, 23)), side = "left")

save.but <- tkbutton(tfr, text = "Save file", width=18, command = function() {
	if(viewer == "TIC" & length(NAME)>0){
		dir.create(paste(tclvalue(textoutput), "/TIC Plot", sep = ""), showWarnings = F)
		if(length(NAME)==1){
			png(paste(tclvalue(textoutput), "/TIC Plot/", NAME, "_RT", as.numeric(tclvalue(RT_begin)), "(", as.numeric(tclvalue(RT_range)), ").png", sep = ""), width = 1440, height = 1080)
		}else{
			png(paste(tclvalue(textoutput), "/TIC Plot/RT", as.numeric(tclvalue(RT_begin)), "(", as.numeric(tclvalue(RT_range)), ")_", Sys.Date(), "-", format(Sys.time(), "%H%M"), ".png", sep = ""), width = 1440, height = 1080)
		}
		TICplot(meta_NAME = NAME, meta_mass = mass, meta_peak = peak, meta_RT = RT, meta_TIC = TIC, rt_begin = as.numeric(tclvalue(RT_begin)), rt_range = as.numeric(tclvalue(RT_range)), cex = 1, inset = c(-0.06, 0.08), cex.axis = 1.5, cex.main = 2, cex.lab = 1.5)
	}else if(viewer == "2D" & length(NAME)>0){
		dir.create(paste(tclvalue(textoutput), "/2D Plot", sep = ""), showWarnings = F)
		png(paste(tclvalue(textoutput), "/2D Plot/", NAME, "_RT", as.numeric(tclvalue(RT_begin)), "(", as.numeric(tclvalue(RT_range)), ")_mz", as.numeric(tclvalue(MZ_begin)), "(", as.numeric(tclvalue(MZ_range)), ").png", sep = ""), width = 1440, height = 1080)
		Meta_2D.plot(meta_mass = mass, meta_peak = peak, meta_RT = RT, rt_begin = as.numeric(tclvalue(RT_begin)), rt_range = as.numeric(tclvalue(RT_range)), m_reso = as.numeric(tclvalue(m_r)), t_reso = as.numeric(tclvalue(t_r)), u_inten = as.numeric(tclvalue(upperbound)), l_inten = as.numeric(tclvalue(lowerbound)), mz_begin = as.numeric(tclvalue(MZ_begin)), mz_range = as.numeric(tclvalue(MZ_range)))
	}else if(viewer == "3D" & length(NAME)>0){
		dir.create(paste(tclvalue(textoutput), "/3D Plot", sep = ""), showWarnings = F)
		png(paste(tclvalue(textoutput), "/3D Plot/", NAME, "_RT", as.numeric(tclvalue(RT_begin)), "(", as.numeric(tclvalue(RT_range)), ")_mz", as.numeric(tclvalue(MZ_begin)), "(", as.numeric(tclvalue(MZ_range)), ")_h", as.numeric(tclvalue(h_an)), "_v", as.numeric(tclvalue(v_an)), ".png", sep = ""), width = 1440, height = 1080)
		Meta_3D.plot(meta_mass = mass, meta_peak = peak, meta_RT = RT, rt_begin = as.numeric(tclvalue(RT_begin)), rt_range = as.numeric(tclvalue(RT_range)), m_reso = as.numeric(tclvalue(m_r)), t_reso = as.numeric(tclvalue(t_r)), u_inten = as.numeric(tclvalue(upperbound)), l_inten = as.numeric(tclvalue(lowerbound)), h_angle = as.numeric(tclvalue(h_an)), v_angle = as.numeric(tclvalue(v_an)), mz_begin = as.numeric(tclvalue(MZ_begin)), mz_range = as.numeric(tclvalue(MZ_range)))
	}
	dev.off()
})
tkpack(save.but, side = "left", anchor = "s")
tkconfigure(save.but, state = "disable")

tkpack(tklabel(tfr, text = "      ", width = switch(ws, 27, 25, 23)), side = "left")
tkpack(tkbutton(tfr, text = "Exit",width=18, command = function() tkdestroy(tt)), side = "left", anchor = "s")

hsc <<- tclVar()
tclvalue(hsc) <- 2
vsc <<- tclVar()
tclvalue(vsc) <- 1.75
tkpack(tfr <- tkframe(plotframe), side = "bottom", fill = "x")
tkpack(tklabel(tfr, text = "", height = 0, font = fontIntro_para), side = "bottom")
tkpack(tklabel(tfr, text = "", height = 0, font = fontIntro_para), side = "top")
tkpack(tklabel(tfr, text = "Hscale: "), side = "left")
tkpack(tkentry(tfr, textvariable = hsc, width = 6), side = "left")
tkpack(tklabel(tfr, text = "      Vscale: "), side = "left")
tkpack(tkentry(tfr, textvariable = vsc, width = 6), side = "left")












    
tkpack(fr <- tkframe(plotframe, relief = "ridge", borderwidth = 3), side = "bottom")



fr1 <- tkframe(fr, relief = "ridge", borderwidth = 3)
tkpack(fr1, side = "left")
tkpack(tklabel(fr1, text = "Resolution"), side = "top", anchor = "nw")

tkpack(fr1.1 <- tkframe(fr1), side = "top")
tkpack(tklabel(fr1.1, text = "m/z resolution"), side = "left", anchor = "s", pady = 10)
m_r <- tclVar()
tclvalue(m_r) <- m_reso
tkpack(tl_m <- tkentry(fr1.1, width = "7", textvariable = m_r), side = "left")

tkconfigure(tl_m, state = "disable")

tkpack(fr1.2 <- tkframe(fr1), side = "top")
tkpack(tklabel(fr1.2, text = "RT resolution  "), side = "left", anchor = "s", pady = 10)
t_r <- tclVar()
tclvalue(t_r) <- t_reso
tkpack(tl_t <- tkentry(fr1.2, width = "7", textvariable = t_r), side = "left")

tkconfigure(tl_t, state = "disable")

tkbind(tl_m, "<Return>", function(...) tkrreplot(img, hscale = as.numeric(tclvalue(hsc)), vscale = as.numeric(tclvalue(vsc))))
tkbind(tl_t, "<Return>", function(...) tkrreplot(img, hscale = as.numeric(tclvalue(hsc)), vscale = as.numeric(tclvalue(vsc))))





fr2 <- tkframe(fr, relief = "ridge", borderwidth = 3)
tkpack(fr2, side = "left")
tkpack(tklabel(fr2, text = "Intensity limit"), side = "top", anchor = "nw")


tkpack(fr2.1 <- tkframe(fr2), side = "top")
tkpack(tklabel(fr2.1, text = "Lower bound  "), side = "left", anchor = "s", pady = 10)
lowerbound <- tclVar()
tclvalue(lowerbound) <- lb_inten
tkpack(tl_lb <- tkentry(fr2.1, width = "7", textvariable = lowerbound), side = "left")


tkpack(fr2.2 <- tkframe(fr2), side = "top")
tkpack(tklabel(fr2.2, text = "Upper bound  "), side = "left", anchor = "s", pady = 10)
upperbound <- tclVar()
tclvalue(upperbound) <- ub_inten
tkpack(tl_ub <- tkentry(fr2.2, width = "7", textvariable = upperbound), side = "left")

tkconfigure(tl_lb, state = "disable")
tkconfigure(tl_ub, state = "disable")
	

tkbind(tl_lb, "<Return>", function(...) tkrreplot(img, hscale = as.numeric(tclvalue(hsc)), vscale = as.numeric(tclvalue(vsc))))
tkbind(tl_ub, "<Return>", function(...) tkrreplot(img, hscale = as.numeric(tclvalue(hsc)), vscale = as.numeric(tclvalue(vsc))))





fr3 <- tkframe(fr, relief = "ridge", borderwidth = 3)
tkpack(fr3, side = "left")
tkpack(tklabel(fr3, text = "Retention time"), side = "top", anchor = "nw")

tkpack(fr3.1 <- tkframe(fr3), side = "top")
tkpack(tklabel(fr3.1, text = "Begin"), side = "left", anchor = "s", pady = 10)
RT_begin <- tclVar()
tclvalue(RT_begin) <- rt_begin

tkpack(ts_RT_begin <- tkscale(fr3.1, width = "14", variable = RT_begin, orient = "horizontal"), side = "left")

ts_RT_begin_23 <- tkscale(fr3.1, width = "14", variable = RT_begin, orient = "horizontal")


tkconfigure(ts_RT_begin, state = "disable")

tkpack(fr3.2 <- tkframe(fr3), side = "top")
tkpack(tklabel(fr3.2, text = "Range"), side = "left", anchor = "s", pady = 10)
RT_range <- tclVar()
tclvalue(RT_range) <- rt_range

tkpack(ts_RT_range <- tkscale(fr3.2, width = "14", variable = RT_range, orient = "horizontal"), side = "left")

ts_RT_range_23 <- tkscale(fr3.2, width = "14", variable = RT_range, orient = "horizontal")

tkconfigure(ts_RT_range, state = "disable")





fr4<-tkframe(fr, relief = "ridge", borderwidth = 3)
tkpack(fr4, side = "left")
tkpack(tklabel(fr4, text = "M/z"), side = "top", anchor = "nw")

tkpack(fr4.1<-tkframe(fr4), side = "top")
tkpack(tklabel(fr4.1, text = "Begin"), side = "left", anchor = "s", pady = 10)
MZ_begin <- tclVar()
tclvalue(MZ_begin) <- mz_begin
tkpack(ts_MZ_begin <-tkscale(fr4.1, width = "14", variable = MZ_begin, orient = "horizontal"), side = "left")

tkconfigure(ts_MZ_begin, state = "disable")

tkpack(fr4.2 <- tkframe(fr4), side = "top")
tkpack(tklabel(fr4.2, text = "Range"), side = "left", anchor = "s", pady = 10)
MZ_range <-tclVar()
tclvalue(MZ_range) <- mz_range
tkpack(ts_MZ_range <- tkscale(fr4.2, width = "14", variable = MZ_range, orient = "horizontal"), side = "left")

tkconfigure(ts_MZ_range, state = "disable")






fr5 <- tkframe(fr, relief = "ridge", borderwidth = 3)
tkpack(fr5, side = "left")
tkpack(tklabel(fr5, text = "Rotation angle"), side = "top", anchor = "nw")

tkpack(fr5.1<-tkframe(fr5), side = "top")
tkpack(tklabel(fr5.1, text = "V angle"), side = "left", anchor = "s", pady = 10)
v_an <- tclVar()
tclvalue(v_an) <- v_angle
tkpack(ts_v_angle <-tkscale(fr5.1, width = "14", variable = v_an, orient = "horizontal"), side = "left")

tkconfigure(ts_v_angle, state = "disable")

tkpack(fr5.2 <- tkframe(fr5), side = "top")
tkpack(tklabel(fr5.2, text = "H angle"), side = "left", anchor = "s", pady = 10)
h_an <-tclVar()
tclvalue(h_an) <- h_angle
tkpack(ts_h_angle <- tkscale(fr5.2, width = "14", variable = h_an, orient = "horizontal"), side = "left")

tkconfigure(ts_h_angle, state = "disable")
tkpack.configure(fr1, fill = "y")

tkpack(frameforplot <- tkframe(plotframe, width = switch(ws, 800, 760, 720), height = switch(ws, 600, 480, 420)), fill = "y", side = "top")


tkfocus(tt)
tkwait.window(tt)
