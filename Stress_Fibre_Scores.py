# Tiled analysis of actin stress fibres. Uses the Kirsch directional derivative estimators to detect the
# principal gradient direction and scores for the presence of stress fibres by calculating the relative contribution
# of this principal direction to the overall intensity gradient.
#
#	-by Richard Butler, Gurdon Institute Imagining Facility


import math as maths

from ij import IJ, WindowManager, ImagePlus, ImageStack, Prefs
from ij.process import ImageProcessor, Blitter
from ij.measure import Calibration, ResultsTable
from ij.gui import Roi, Arrow, TextRoi, Overlay

from java.awt import Color, Font
from java.io import File

from loci.plugins import BF
from loci.plugins.in import ImporterOptions


font = Font(Font.SANS_SERIF, Font.BOLD, 12)
tileRoiColour = Color(1.0, 1.0, 0.0, 0.25)

#Kirsch kernels
dirKernels = [
	[5,5,5,
	-3,0,-3,
	-3,-3,-3],
	[-3,5,5,
	-3,0,5,
	-3,-3,-3],
	[-3,-3,5,
	-3,0,5,
	-3,-3,5],
	[-3,-3,-3,
	-3,0,5,
	-3,5,5],
	[-3,-3,-3,
	5,0,-3,
	5,5,-3],
	[-3,-3,-3,
	-3,0,-3,
	5,5,5],
	[5,5,-3,
	5,0,-3,
	-3,-3,-3],
	[5,-3,-3,
	5,0,-3,
	5,-3,-3]
]

def process(inip):
	ip = inip.duplicate()
	sub = inip.duplicate()
	sigmaPx = 2.0
	ip.blurGaussian(sigmaPx)
	sub.blurGaussian(sigmaPx*1.4)
	ip.copyBits(sub, 0,0, Blitter.SUBTRACT)
	return ip

def normAngle(theta):
	if theta<0:
		theta = 2*maths.pi + theta
	if theta > maths.pi:
		theta = theta - maths.pi
	return theta

def edgeDirections(ip):
	values = []
	for k in dirKernels:
		d = ip.duplicate()
		d.convolve(k, 3, 3)
		values.append(d.getStatistics().mean)
	return values

def tiledAnalysis(ip, ol, title):
	W = ip.getWidth()
	H = ip.getHeight()
	if W!=H:
		IJ.error("Image is not square")
		return None

	nRows = 24	#number of rows and columns of tiles

	row = rt.getCounter()
	scores = []
	step = int(H/nRows)	#px
	overlap = int(step/5.0)
	for y in range(overlap,H-1,step):
		for x in range(overlap,W-1,step):
			ox = x+step/2
			oy = y+step/2
			if ox>W-1 or oy>H-1: continue
			
			tileRoi = Roi(x-overlap,y-overlap, step+overlap*2, step+overlap*2)
			ip.setRoi(tileRoi)
			tile = ip.crop()
			tileStats = tile.getStatistics()
			fibreI = tile.getStatistics().mean

			edgeDir = edgeDirections(tile)
			edgeSum = sum(edgeDir)
			edgeMax = max(edgeDir)
			edgeMin = min(edgeDir)
			score = ((edgeMax-edgeMin)/(edgeSum-edgeMin)) * 8
			maxi = edgeDir.index(edgeMax)
			theta = maxi*((2*maths.pi)/float(len(edgeDir))) + (maths.pi/2.0)
			theta = normAngle(theta)

			rt.setValue("Image", row, title)
			rt.setValue("Tile", row, row)
			rt.setValue("X", row, ox*cal.pixelWidth)
			rt.setValue("Y", row, oy*cal.pixelHeight)
			rt.setValue("Fibre Intensity", row, fibreI)
			rt.setValue("Principal Direction", row, theta)
			rt.setValue("Stress Fibre Score", row, score)
			
			tileRoi.setStrokeColor(tileRoiColour)
			ol.add(tileRoi)
			
			scale = 20	#scaling factor for arrow display
			dx = maths.cos(theta)*scale
			dy = maths.sin(theta)*scale

			arrow = Arrow(ox-dx, oy+dy, ox+dx, oy-dy)
			arrow.setStyle(Arrow.OPEN)
			arrow.setDoubleHeaded(True)
			arrow.setStrokeWidth(4)
			arrow.setHeadSize(4)
			f = max(0.0, min(score, 1.0))
			arrowColour = Color(1-f, f, f/2.0)
			arrow.setStrokeColor(arrowColour)
			ol.add(arrow)

			label = TextRoi(x+step/2-10, y+6, str(row), font)
			label.setStrokeColor(Color.CYAN)
			ol.add(label)

			row += 1

def run(imp):
	global cal
	cal = imp.getCalibration()
	stack = imp.getImageStack()
	actin = None
	maxM = -1
	maxZ = -1
	for z in range(imp.getNSlices()):
		ip = stack.getProcessor( imp.getStackIndex(2, z, 1) )
		value = ip.getStatistics().mean
		if value > maxM:
			maxM = value
			maxZ = z
	actin = stack.getProcessor( imp.getStackIndex(2, maxZ, 1) )
	
	ol = Overlay()
	title = imp.getTitle()+"-Z"+str(maxZ)
	proc = process(actin)
	
	tiledAnalysis(proc, ol, title)

	actinImp = ImagePlus(title, actin)
	actinImp.setCalibration(cal)
	IJ.run(actinImp, "Grays", "")
	actinImp.setOverlay(ol)
	actinImp.show()

def fileHandler(fil):
	if fil.isFile():
		path = fil.getAbsolutePath()
		if path.endswith(".lif"):
			options = ImporterOptions()
			options.setId(path)
			options.setOpenAllSeries(True)
			images = BF.openImagePlus(options)
			for imp in images:
				run(imp)
				imp.close()
	elif fil.isDirectory():
		files = fil.listFiles()
		for f in files:
			fileHandler(f)

rt = ResultsTable()
rt.showRowNumbers(False)
if WindowManager.getCurrentImage() == None:
	path = IJ.getDirectory("lif or Directory")
	fileHandler(File(path))
	rt.show("Stress Fibres : "+path)
else:
	imp = IJ.getImage()
	run(imp)
	rt.show("Stress Fibres : "+imp.getTitle())
