
#this script is used to split the process of Ursa Major with some fundamental operations of FITS

from myPYTHON import *
#import table

from astropy.table import Table
#print dir(plt)
from matplotlib.patches import Polygon
class UrsaMajor:

	"""
	This class is used to deal with all UrsaMajor issues
	"""
	#initialize

	originalFITSPath="./OriginalFITS/"

	FITS12COName="myCut12CO.fits" 
	FITS13COName="myCut13CO.fits" 
	FITSC18OName="myCutC18O.fits" 

	FITS12CO=originalFITSPath+"myCut12CO.fits" 
	FITS13CO=originalFITSPath+"myCut13CO.fits" 
	FITSC18O=originalFITSPath+"myCutC18O.fits" 


	middleOutPutPath="./middleFITS/"

	tempOutPutPath="./tempFITS"

	data12CO,header12CO=None,None
	data13CO,header13CO=None,None

	dataC18O,headerC18O=None,None


	doFITS= None #myFITS()

	HIFITS=originalFITSPath+"HICube.fits"
	IRASFITS=originalFITSPath+"IRAS100.fits"

	distance=110 #pc
	Vrange=[-6,9] #mark the range for moment km/s



	def __init__(self,read12CO=True,read13CO=True,readC18O=False):

		#read 12CO and 13CO
		#no siginals from 18O, we will see, examine the peak positions of C18O to examine its limitition
				
		self.doFITS=myFITS()
		if read12CO:
			data12CO,header12CO=self.doFITS.readFITS(self.FITS12CO)
		if read13CO:
			data13CO,header13CO=self.doFITS.readFITS(self.FITS13CO)
		if readC18O:
			dataC18O,headerC18O=self.doFITS.readFITS(self.FITSC18O)
		#myFITS.getSurvey("IRIS 100",LB=[144.92432,38.395361],Pixels=800)
		#myFITS.getSurvey("IRIS 60",LB=[144.92432,38.395361],Pixels=800)

	def draw12CORGB(self,regrid=False):
		self.drawRGT(self.data12CO,self.header12CO,regrid=regrid)

	def draw13CORGB(self,regrid=False):
		self.drawRGT(self.data13CO,self.header13CO,regrid=regrid)

	def drawC18ORGB(self,regrid=False):
		self.drawRGT(self.dataC18O,self.headerC18O,regrid=regrid)





	def addShadow(self,s1,ax1,redVs,colorCode):
		"""
		What the hell is zeroPosition
		"""
		xdata=s1.get_xdata()
		ydata=s1.get_ydata()
		
		a=float(redVs[0]);b=float(redVs[1])
		
		indexa=0;indexb=0
		
		for x in xdata:
			if x<a:
				indexa=indexa+1
			if x<b:
				indexb=indexb+1
			if x>b:
				break
				
		iy=ydata[indexa:indexb+1]
		ix=xdata[indexa:indexb+1]
			
		verts=[(xdata[indexa],0)]+list(zip(ix,iy))+[(xdata[indexb],0)]
		
		
		
		#print verts
		poly=Polygon(verts,facecolor=colorCode,edgecolor='0.75',linewidth=0.2)
		ax1.add_patch(poly)
		



	def getSpectraByIndex(self,data,dataHeader,indexX,indexY):

		"""
		This function is used to get a voxel value from a data cube
		
		v, km/s
		"""
		wcs = WCS(dataHeader)
 
 
		
		##it is possible that the yindex,and xindex can exceed the boundary of spectrum
		
 
 
		spectral=data[:, indexY,indexX]
		##just for test
		#print  w.all_world2pix(0, 0, 0,0)
		#print data.shape[0]
		velocityIndex= range(data.shape[0])
		
		velocities=wcs.all_pix2world(0, 0,velocityIndex,0)[2]/1000.
 
		# 
		return spectral,velocities
	def getAverageSpec(self,fitsFile,path="./"):
		"""
		This function is dedicated to get the average spectrum for the CO lines
		
		#Basicly, we only care about the spectra at the peak posion of the cores
		
		#The index of peak position given by Duchamp is from 0.
		
		"""
		#read the file
		cores=self.getCoreTable()

		COdata,COheader=self.doFITS.readFITS(path+fitsFile)

		#print len(cores)
		avgSpec=0
 
		for eachCore in cores:
			#l,b= eachCore["GLON_deg"],eachCore["GLAT_deg"]
			#spectrum,vs =self.getSpectraByLB( COdata,COheader,l,b) 
			
			X,Y= eachCore["X_peak"],eachCore["Y_peak"]
			spectrum,vs =self.getSpectraByIndex( COdata,COheader,int(X),int(Y)) 
			
			
			avgSpec=avgSpec+spectrum
			#print l,b,spectrum[0]
		avgSpec=avgSpec/1./len(cores)
 
		if 0:
			l,b= cores[0]["GLON_deg"],cores[0]["GLAT_deg"]
	
			avgSpec,vs =self.getSpectraByLB( COdata,COheader,l,b) 
		
		if 0:
			fig, ax = plt.subplots()
			ax.plot(vs,avgSpec)
			plt.show()

		return avgSpec,vs


	def getCoreTable(self):
		
		"""
		This function get the table of cores, and eliminate the edge fake cores
		
		"""


		duchampFile="/Users/qzyan/WORK/NewUrsaMajorPaper/Duchamp/duchamp-Results.txt"
		
		colNamesStr=" ObjID            Name      X      Y     Z       GLON      GLAT       GLON_deg       GLAT_deg      VRAD      MAJ      MIN     PA   w_GLON   w_GLAT     w_50     w_20   w_VRAD        F_int   F_peak   X2  Y1  Y2  Z1  Z2 Nvoxel Nchan Nspatpix Flag   X_av   Y_av  Z_av X_cent Y_cent Z_cent X_peak Y_peak Z_peak"
		
		from astropy.io import ascii
		
		tempFile="/Users/qzyan/WORK/NewUrsaMajorPaper/Duchamp/duchamp-Results.dat"
		##This step is to filter the lines Commented.
		#
		
		#The file used to identify cores
		#CO13FITS="/home/qzyan/HLM/FITS/Duchamp/finalALL_C.fits"
		
		CO13data,CO13Header=self.doFITS.readFITS(self.FITS13CO)
		#
		#header=fits.getheader(CO13FITS)
		#
		#WCS(header)
 
		
 
		
		#print WCS(CO13Header)
		
		data_rows=[]
		
		
		for line in open(duchampFile):
			if line.startswith('#'):
				
			  continue
			 
			if line.strip()=="":
				continue
 
			data_rows.append(line.split())
 
		
		colNames=colNamesStr.split()
		
 
		
		dtypes=[]
		#print 
		for eachValue in data_rows[0]:
			
			try :
				float(eachValue)
				dtypes.append("f8") 
			except:
				dtypes.append('S50')
			
			
			#print eachValue, eachValue.replace('.','',1).isdigit()
		
 
		rawCores=Table(rows=data_rows,names=colNames,dtype=dtypes)
 
		
		
		myCores=rawCores[0:1]
	
		myCores.remove_row(0)
		
		indexCore=1
		
		for eachCore in rawCores:
			
			#filter the cores, eliminate the fake cores (edge fake)
			l,b=eachCore["GLON_deg"],eachCore[ "GLAT_deg"]
			dL=2.5/60.
			#check if with in 2 arcmin, nan value exist
			value1=self.doFITS.getVoxValue(CO13data,CO13Header,[l+dL,b,0])
			value2=self.doFITS.getVoxValue(CO13data,CO13Header,[l-dL,b,0])
			value3=self.doFITS.getVoxValue(CO13data,CO13Header,[l,b-dL,0])
			value4=self.doFITS.getVoxValue(CO13data,CO13Header,[l,b+dL,0])
			
 
			
			if  np.isnan(value1) or np.isnan(value2)  or np.isnan(value3)  or np.isnan(value2):
				#print indexCore
				
				indexCore+=1
			#print value1,value2,value3,value4
				continue
			
			#elset 
			myCores.add_row(eachCore)
			indexCore+=1

		return myCores





	def drawCOInt(self,fitsFile):
		
		"""
		This fucntion is dedicatted to to plot of RGB image with different velocity range as different colors
		
		"""
	
		#first test one color
	
		#range1=[-6,8] # /km/s # negative blue, positive red

		range1=[-7,-1]
		range2=[-1,2]
		range3=[2,9]
			

 

			
		#f = pyfits.open(fitsFile)
		#c = cube.Cube(f[0].data, f[0].header)
		
		
		bData,header1=self.doFITS.momentFITS(fitsFile,range1,0,"./tempFITS/COFalse.fits")
 		if bData.shape[0]==1:
			bData=bData[0]
 
		gData,header2=self.doFITS.momentFITS(fitsFile,range2,0,"./tempFITS/COFalse.fits")
 		if gData.shape[0]==1:
			gData=gData[0]
	 
		
		rData,header3=self.doFITS.momentFITS(fitsFile,range3,0,"./tempFITS/COFalse.fits")
 		if rData.shape[0]==1:
			rData=rData[0]

 
 
		
		img = np.zeros((bData.shape[0], bData.shape[1], 3), dtype=float)

		scales=[0.1,8]

		if "13" in fitsFile:
			scales=[0.01,2]


		img[:,:,0]=img_scale.sqrt(rData, scale_min=scales[0],scale_max=scales[1])
		img[:,:,1]= img_scale.sqrt(gData,scale_min=scales[0],scale_max=scales[1])
		img[:,:,2]=img_scale.sqrt(bData, scale_min=scales[0],scale_max=scales[1])



		drawWcs=WCS(header1)
		
		ax=pywcsgrid2.subplot(111,header=drawWcs)
		#ax.imshow(bData, origin="lower")
		#from astropy.visualization import make_lupton_rgb
		#image = make_lupton_rgb(rData, gData, bData, stretch=0.5)
		rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
		rc('text', usetex=True)
		#ax.imshow(bData,origin="lower",cmap="jet")
		i = imshow_rgb(ax, self.doFITS.mytrim(rData,0,8) ,self.doFITS.mytrim(gData,0,8) , self.doFITS.mytrim(bData,0,8) , origin="lower", interpolation="none",vmin=2,vmax=15 )

		#i = imshow_rgb(ax, img[:,:,0] ,img[:,:,1] ,img[:,:,2] , origin="lower", interpolation="none",vmin=2,vmax=15 )
		#axins = zoomed_inset_axes(ax,   loc=3)
		axins=inset_axes(ax, width="32%", height="25%", loc=1)
		
		#axins=plt.axes([0.3,0.17,0.2 ,0.18]) 
		axins.set_axis_bgcolor('lightcyan')
		
		spectral12,vs12=self.getAverageSpec(self.FITS12COName,path=self.originalFITSPath)
		spectral13,vs13=self.getAverageSpec(self.FITS13COName,path=self.originalFITSPath)



		s1,=axins.step(vs12,spectral12,color='black',linewidth=0.3,label=r"$^{12}$CO")
		##add shadow
		self.addShadow( s1,axins,range1,"blue")
		self.addShadow( s1,axins,range2,"green")
		self.addShadow( s1,axins,range3,"red")
		
		
		axins.step(vs13,spectral13,color='purple',linewidth=0.3,label=r"$^{13}$CO")
		
		axins.set_ylim([-0.2,1.6])
		axins.set_xlim([min(vs12),max(vs12)])
		#axins.tick_params(axis='both', which='major', labelsize=5,direction='in',colors="white", pad=-12)
		axins.tick_params(axis='both', which='major', labelsize=5,direction='out',colors="white",top=False,right=False )

		axins.set_xlabel(r"Velocity (km s$^{-1}$)", fontsize=5 )
		axins.set_ylabel("Brightness temperature (K)", fontsize=5)

		
		
		axins.xaxis.label.set_color('white')
		axins.yaxis.label.set_color('white')

		legend=axins.legend(fontsize=5)
		legend.get_frame().set_linewidth(0.2)
		## draw base line
		legend.get_frame().set_facecolor('paleturquoise')
		axins.plot(vs12,vs12*0,color='black',linewidth=0.2)

		#axins.set_label_position('top')
		#axins.xaxis.label.set_color('white')
		#axins.yaxis.label.set_color('white')
		#axins.tick_params(direction='in')
		#axins.set_axislabel_direction("+")
		
		
		# HI, pix resolution 
		HIPix=0.5 # arcmin
		resoltion=np.radians(HIPix/60)
 
		#what is the angle of 1 pc
 
		angle1PC=2*np.arctan(0.5/self.distance)
		size1PC=angle1PC/resoltion

  
		ax.add_size_bar(size1PC, # 30' in in pixel
								r"1 pc",
								loc=8,color='paleturquoise')

		
		#.tick_params(axis='both', which='minor', labelsize=8)
		
		#plt.show()
		#ax.imshow( )
		#calculate the spectral
		
		ax.axis[:].major_ticks.set_color("w")
		
		#ax.locator_params(axis="x", nbins=6)
		#ax.add_inner_title("Figure 1", loc=2)
		#plt.show()
		

		if "12" in fitsFile:
			saveName="RGBInt12CO.pdf"
		if "13" in fitsFile:
			saveName="RGBInt13CO.pdf"
		if "18" in fitsFile:
			saveName="RGBIntC18O.pdf"

		plt.savefig(saveName, bbox_inches="tight")
		#plt.show()
		plt.clf()
		#
	




	def drawRGB(self,FITSdata,FITSheader,regrid=False):
		"""
		This function is dedicated to draw RGB of CO, HI, and iras 100 um.
		
		#regrided by miriad
		
		#/home/qzyan/HLM/FITS/project/CORGB.fits
		#/home/qzyan/HLM/FITS/project/HIRGB.fits
		#/home/qzyan/HLM/FITS/project/IRASRGB.fits
		
		
		/home/qzyan/HLM/rgd12CO.fits
		/home/qzyan/HLM/rgdIRAS100.fits
		/home/qzyan/HLM/FITS/UrsaMjaorHI.fits
		#Co and IRAS have been smoothed and resemble, Tog get the RGB fits, moment is needted
		
		"""
		
		# moment fits
		
		COInt="CORGB.fits"
		HIInt="HIRGB.fits"
		
		#"./FITS/CO12.fits"
		
		CO12Data,CO12Header=self.readFITS(CO12FITS)
		#do smooth first
		if regrid:
			self.smoothSpaceFITS(CO12Data, CO12Header,52/60.,16.2,"Smooth12CO.fits",) # the resolution of HI is 16.2 arcmin
			
			self.regridFITS( "Smooth12CO.fits","./FITS/UrsaMjaorHI.fits",'rgd12CO.fits')
		
		
		self.momentFITS("rgd12CO.fits",[-6,9],0,"./",COInt)
		self.momentFITS("UrsaMjaorHI.fits",[-6,9],0,"./FITS/",HIInt)
		
		
		
		rData,rHeader=self.readFITS("./FITS/"+HIInt)
		gData,gHeader=self.readFITS("./"+COInt)
		bData,bHeader=self.readFITS("./rgdIRAS100.fits")
		
		print rData.shape
		print gData.shape
		print bData.shape
		
 
		
		img = np.zeros((rData.shape[1], rData.shape[2], 3), dtype=float)
		
		
 
		
		img[:,:,0] = img_scale.sqrt(rData[0], scale_min=120, scale_max=266)
		img[:,:,1] = img_scale.sqrt(gData[0], scale_min=0.5, scale_max=6)
		img[:,:,2] = img_scale.sqrt(bData , scale_min=3, scale_max=5.5)
		
		ax=pywcsgrid2.subplot(111,header=rHeader)
		#ax.imshow(bData, origin="lower")
		#$rc('text', usetex=True)
		#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
		rc('font',**{'family':'serif','serif':['Times']})
 		i = imshow_rgb(ax, img[:,:,0] ,img[:,:,1] , img[:,:,2] , origin="lower", interpolation="none",vmin=2,vmax=15 )
		ax.axis[:].major_ticks.set_color("w")


		pgf_with_latex = {
		    "text.usetex": True,            # use LaTeX to write all text
		    "pgf.rcfonts": False,           # Ignore Matplotlibrc
		    "pgf.preamble": [
		        r'\usepackage{color}'     # xcolor for colours
		    ]
		}
		#pl.rcParams.update(pgf_with_latex)
		rc('text', usetex=True)
		#matplotlib.use('ps')
		fp = dict(size=12  )
		rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
		rc('text', usetex=True)
		rc('text.latex', preamble='\usepackage{color}')
		
		
		import matplotlib
		from matplotlib.backends.backend_pgf import FigureCanvasPgf
		matplotlib.backend_bases.register_backend('pdf', FigureCanvasPgf)
		plt.rcParams.update(pgf_with_latex)
		
		#plt.rcParams["text.usetex"] = True #mpl.rcParams['text.usetex'] = True
		#cName = AnchoredText(r"(R, G, B)= (\textcolor{blue}{ H\,\textsc{i}},$^{12}$CO ($J=1\rightarrow 0$), IRAS 100 $\mu$m)",loc=3, prop=fp,frameon=False,  )

		cName = AnchoredText(r"\textcolor{red}{H\,\textsc{i} }\textcolor{green}{$^{12}$CO ($J=1\rightarrow 0$) } \textcolor{blue}{IRAS 100 $\mu$m}",loc=3, prop=fp,frameon=False,  )

		ax.add_artist(cName)		

		#plt.show()
		ax.set_xlim(140,258)
		ax.set_ylim(56,171)

 
		# HI, pix resolution 
		HIPix=5 # arcmin
		resoltion=np.radians(5./60)
 
		#what is the angle of 1 pc
 
		angle1PC=2*np.arctan(0.5/self.distance)
		size1PC=angle1PC/resoltion
 		
		ax.add_size_bar(size1PC, # 30' in in pixel
								r"\textcolor{white}{1 pc}",
								loc=1,color='w')
		

		#pywcsgrid2.savefig('RGBCIHIIRAS.eps', bbox_inches="tight")
		plt.savefig('RGBCOHIIRAS.pdf', bbox_inches="tight")

		#print "???"
		
	



	def printStatus(self):
		pass

    





doHL=UrsaMajor()


if 0: #Examine a possibile outflow
	#doHL.examineOutflow()
	pass

	#downLoadWiseData
	#doHL.doFITS.detectSurveyResolution("WISE 22",[148,38])
	#doHL.doFITS.getSurvey("WISE 22",[148,38],Pixels=1000)
	#doHL.doFITS.downLoadSurveyByRange("WISE 22",LRange=[138,152],BRange=[30,46],size=1.,Pixels=1000)


if 0: #drawRGB of CO IRAS HI
	doHL.drawRGB("./FITS/CO12_C.fits",regrid=False)


if 1: #draw 12CO three false color 
	doHL.drawCOInt("./OriginalFITS/myCut12CO.fits")

	doHL.drawCOInt("./OriginalFITS/myCut13CO.fits")

if 0:
	pass
	#doHL.drawRGB()
	#doHL.draw12CO("CO12.fits")
	#doHL.momentFITS("CO12.fits",[-6,9],0,"./FITS/")
	

	doHL.doFITS.momentFITS("./OriginalFITS/myCut12CO.fits",[-6,8.2],0,"./middleFITS/12COint.fits")

	#doHL.momentFITS("UrsaMjaorHI.fits",[-6,8],0,"./FITS/")

 

if 0: # testing the average spectrum
	avgSp12=doHL.doFITS.getAverageSpec("CO13.fits",path="./FITS/")
	fig, ax = plt.subplots()
	
	#print avgSp12
 	ax.plot(avgSp12 )
	plt.show()


if 0: #tesing Duchump, using duchump to find the core of 13CO, 
	pass
	
if  0: # testing resemble
	pass
	#doHL.momentFITS(FITSFile,Vrange,mom,outPUTPATH,cuttEdge=False):##Vrange kms
	
	
	# smooth CO data
	#CO12Data,CO12Header=doHL.readFITS("./FITS/CO12_C.fits")
	
	doHL.smoothSpaceFITS(CO12Data, CO12Header,52/60.,16.2,"Smooth12CO.fits") # the resolution of HI is 16.2 arcmin
	# regrid CO data
	#doHL.regridFITS("Smooth12CO.fits","./FITS/UrsaMjaorHI.fits",'rgd12CO.fits')
	
	
	
	## smooth IRAS data
	#IRAS100Data,IRAS100Header=doHL.readFITS("./FITS/IRAS100.fits")
	#doHL.smoothSpaceFITS(IRAS100Data, IRAS100Header,3,16.2,"SmoothIRAS100.fits") # the resolution of HI is 16.2 arcmin

	
	#doHL.regridFITS("SmoothIRAS100.fits","./FITS/UrsaMjaorHI.fits",'rgdIRAS100.fits')

	
if 0: # testing the Cores, mark the cores on 13CO map, calculate the mass for each core, and draw the CMF.
	
	"""
	Testing if core are gravitionaly bound, and 
	"""
	doHL.createClean13COMoment("./FITS/CO13.fits")

	
	pass
	
	# before identifying 13CO cores, we need to clear the 13CO spectral, remove spectral lines which have large RMS
	#createClean13COMoment


if 0: #check the mass distribution of cores identified by Duchamp (13CO)
		
	"""
	333
	"""
	
	doHL.drawCoreMass()
	pass
	
	

if 0: #check the mass distribution of cores identified by Duchamp (13CO)
	
	"""
	333
	"""
	
	doHL.getAverageSpec( "CO13.fits",path="./FITS/") 
	pass
	
	
if 0: #check the RMS of 13CO
	pass
	sigma13=doHL.getRMS("./FITS/CO13.fits")
	print "the RMS of 13CO is ",sigma13
	
	
if 0:
	pass
	#
	doHL.compare12And13()





