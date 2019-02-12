# -*- coding: utf-8 -*-
from __future__ import division

import math
import xc_base
import geom
import xc
from materials import typical_materials as tm
from materials.prestressing import prestressed_concrete as presconc
from rough_calculations import ng_prestressed_concrete as ng_presconc
from model import predefined_spaces
from solution import predefined_solutions
from model.mesh import finit_el_model as fem
from actions import loads
from materials.ehe import EHE_materials

# Units: m, N

#Data
#Geometry
wWeb=0.060
hWeb=0.18
wFlange=0.03
hFlange=0.036+0.018/2
areaTopTendons=2*math.pi*(2.25e-3)**2/4
areaBotTendons=3*math.pi*(2.25e-3)**2/4
Lbeam=5.75
esize=5.75/20
#Materials
fptk=1800e6
Ep=1.95e11
sigma_pi=1200e6
#Loads
g1=330
g2=500

# ****
# Problem type
FEcase= xc.FEProblem()
prep=  FEcase.getPreprocessor
points=prep.getMultiBlockTopology.getPoints
lines=prep.getMultiBlockTopology.getLines
sets=prep.getSets 
nodes= prep.getNodeHandler
elems=prep.getElementHandler
modelSpace= predefined_spaces.StructuralMechanics3D(nodes)

#Materials
concr=EHE_materials.HP45
concrDiag=concr.defDiagK(prep)
dgDHP45= EHE_materials.HP45.getDiagK(prep)

steelP=EHE_materials.EHEPrestressingSteel(steelName='steelP',fpk=fptk,alpha= 0.7, steelRelaxationClass=1, tendonClass= 'strand')
steelPDiag=steelP.defDiagK(prep,sigma_pi)
steelPDiag.E=Ep

#Cross-sectional geometry
geomSect= prep.getMaterialHandler.newSectionGeometry('geomSect')
regions= geomSect.getRegions
#Web
concrSect=regions.newQuadRegion(concrDiag.getName())
concrSect.pMin=geom.Pos2d(-wWeb/2,0)
concrSect.pMax=geom.Pos2d(wWeb/2,hWeb)
concrSect.nDivIJ=8
concrSect.nDivJK=24
#Left flange
concrSect=regions.newQuadRegion(concrDiag.getName())
concrSect.pMin=geom.Pos2d(wWeb/2,0)
concrSect.pMax=geom.Pos2d(wWeb/2+wFlange,hFlange)
concrSect.nDivIJ=3
concrSect.nDivJK=4
#Right flange
concrSect=regions.newQuadRegion(concrDiag.getName())
concrSect.pMin=geom.Pos2d(-wWeb/2-wFlange,0)
concrSect.pMax=geom.Pos2d(-wWeb/2,hFlange)
concrSect.nDivIJ=3
concrSect.nDivJK=4

#Active reinforcement
reinf= geomSect.getReinfLayers
reinfLayer=reinf.newStraightReinfLayer(steelPDiag.getName())
reinfLayer.numReinfBars= 2
reinfLayer.barArea=areaTopTendons
reinfLayer.p1=geom.Pos2d(0,hWeb-0.06)
reinfLayer.p2=geom.Pos2d(0,hWeb-0.02)
reinfLayer=reinf.newStraightReinfLayer(steelPDiag.getName())
reinfLayer.numReinfBars= 1
reinfLayer.barArea=areaBotTendons
reinfLayer.p1=geom.Pos2d(0,0.04)
reinfLayer.p2=geom.Pos2d(0,0.04)

reinfLayer=reinf.newStraightReinfLayer(steelPDiag.getName())
reinfLayer.numReinfBars= 5
reinfLayer.barArea=areaBotTendons
reinfLayer.p1=geom.Pos2d(-wWeb/2-wFlange+0.02,0.02)
reinfLayer.p2=geom.Pos2d(wWeb/2+wFlange-0.02,0.02)
                         
#Section material 
materiales= prep.getMaterialHandler
sctFibers= prep.getMaterialHandler.newMaterial("fiber_section_3d","sctFibers")
fiberSectionRepr= sctFibers.getFiberSectionRepr()
fiberSectionRepr.setGeomNamed("geomSect")
sctFibers.setupFibers()

'''                         
#Plot
from materials.sections.fiber_section import plot_fiber_section as pfs
fsPlot=pfs.fibSectFeaturesToplot(fiberSection=geomSect)
fig1,ax2d=fsPlot.generatePlot()
fig1.show()
fig1.savefig('fig1.png')
'''

#Nodes
l= 1e-7     # Distance between nodes
nodes.defaultTag= 1 #First node number.
nod= nodes.newNodeXYZ(1.0,0,0)
nod= nodes.newNodeXYZ(1.0+l,0,0)

# Element definition
elems.defaultMaterial= "sctFibers"
elems.dimElem= 1 # Dimension of element space
elems.defaultTag= 1 #Tag for the next element.
elem= elems.newElement("ZeroLengthSection",xc.ID([1,2]))
# Constraints
modelSpace.fixNode000_000(1)
modelSpace.fixNodeF00_0F0(2)

sccEl= elem.getSection()         
fibersSccEl= sccEl.getFibers()

#Creation of two separate sets of fibers: concrete and reinforcement steel 
from materials.sections.fiber_section import fiber_sets
setsRCEl= fiber_sets.fiberSectionSetupRCSets(scc=sccEl,concrMatTag=concr.matTagK,concrSetName="concrSetFbEl",reinfMatTag=steelP.matTagK,reinfSetName="reinfSetFbEl")
concrSetFbEl=setsRCEl.concrFibers.fSet
steelSetFbEl=setsRCEl.reinfFibers.fSet
Ap=steelSetFbEl.getArea(1) #area active reinforcement
COGs=steelSetFbEl.getCenterOfMass()
COG=fibersSccEl.getCenterOfMassHomogenizedSection(1)
ecc_p=COG[1]-COGs[1]   #eccentricity
#Pérdida relajamiento t=48h
deltaSgp=steelP.getRelaxationStressLossT(tDays=2, initialStress=sigma_pi)
Sgp48h=sigma_pi-deltaSgp
Np48h=-Sgp48h*Ap
Mp48h=Np48h*ecc_p

for layer in reinf:
    for bar in layer.getReinfBars:
        m=bar.getMaterial()
        m.initialStress+=-deltaSgp

# Self-weight bending moment at mid-span
#Msw=-g1*Lbeam**2/8
        
# Loads definition
loadHandler= prep.getLoadHandler   #loads container
lPatterns= loadHandler.getLoadPatterns
#Load modulation.
ts= lPatterns.newTimeSeries("constant_ts","ts")
lPatterns.currentTimeSeries= "ts"
#Load case definition
lp0= lPatterns.newLoadPattern("default","0")
lPatterns.addToDomain("0") 

# Solve
analisis= predefined_solutions.simple_static_linear(FEcase)
analOk= analisis.analyze(1)
nodes.calculateNodalReactions(True,1e-6)

P_result=steelSetFbEl.getTensionResultant()
P_targ=Ap*(sigma_pi-deltaSgp)
error=(P_targ-P_result)/P_targ
print 'error=',error*100, '%'

'''
# #plot cross-section strains and stresses 
from postprocess import utils_display
# 
#   fiberSet: set of fibers to be represented
#   title:    general title for the graphic
#   fileName: name of the graphic file (defaults to None: no file generated)
#   nContours: number of contours to be generated (defaults to 100)
# 
#utils_display.plotStressStrainFibSet(fiberSet=fibersSccEl,title='cross-section fibers',fileName='problem.jpeg')
utils_display.plotStressStrainFibSet(fiberSet=concrSetFbEl,title='cross-section concrete fibers',fileName='problem.jpeg',nContours=0,pointSize=25,fiberShape='s')
utils_display.plotStressStrainFibSet(fiberSet=steelSetFbEl,title='cross-section steel fibers',fileName='problem.jpeg',nContours=0,pointSize=50)
'''    

