# -*- coding: utf-8 -*-
from __future__ import division
'''Test for checking the deflections in a prestressed concrete beam.
Data for the problem and approximate calculation are taken from 
Example 4.1 of the topic 4 of course "Prestressed Concrete Design 
(SAB 4323) by Baderul Hisham Ahmad 
ocw.utm.my

Problem statement:
Determine the midspan deflection of a beam: 
(i) at transfer with an inertial prestress force of 6800kN; 
(ii) under an imposed load of 30 kN/m when the prestress force has been 
reduced to 4500 kN. 
Take self weight of beam = 11.26 kN/m; I =0.06396m4 ; E = 28 x 10^6 kN/m2
'''

'''Try 3: I set initial stress in prestressing steel =  6800 kN/m2.
Before calculating I try to impose a strain in order to reduce the stress
from 6800 kN/m2 to 4500 kN/m2 (lines 157-162)
results: idem to case 1. The following error arises: ElementBodyLoad::applyLoad; el number of pointers no coincide con el de identifiers.
'''

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

#DATA
#Geometry
span=24      #span of the beam [m]
hBeam=1.305406 #height of the cross-section [m]. Parallel to local z-axis 
wBeam= 0.3450267 #width of the cross-section [m]. Parallel to local y-axis
Abeam=hBeam*wBeam   #cross-section area of the beam[m2]
Iybeam=1/12.*hBeam*wBeam**3 #moment of inertia of the beam cross-section [m4]
Izbeam=1/12.*wBeam*hBeam**3 #moment of inertia of the beam cross-section [m4]
deltaTendon=0.26
nDivLines=8      #number of elements in each line

#Material properties
Ec=28e6      #modulus of elasticity of concrete [kPa]
Ep=Ec*7.5    #modulus of elasticity of prestressing steel
Ep=Ec*1e-3   #modulus of elasticity of prestressing steel (allow big deflection)
nuc=0.2      #coefficient of Poisson of concrete
densc= 2.5   #specific mass of concrete (t/m3)

fy= 1171e3 # Yield stress of the steel expressed in kPa.
Aps=1  #area of tendon cross-section [m2]
#Prestress
fpi=6800/Aps       #initial stress in the tendon [kPa]
fps=4500/Aps       #stress in the tendon in service [kPa]
#Loads
Wsw=Abeam*densc*10   #self-weight uniform load on beam [kN/m]
#END DATA


# XC model of the beam
# Problem type
FEcase= xc.FEProblem()
prep=  FEcase.getPreprocessor
points=prep.getMultiBlockTopology.getPoints
lines=prep.getMultiBlockTopology.getLines
sets=prep.getSets 
nodes= prep.getNodeHandler
modelSpace= predefined_spaces.StructuralMechanics3D(nodes)

#Points and lines beam
beamSet=sets.defSet('beamSet')
beamPoints= beamSet.getPoints
beamLines= beamSet.getLines
for i in range(4):
    p=points.newPntFromPos3d(geom.Pos3d(0,i*span/3,0))
    beamPoints.append(p)
    if(i>0):
        l= lines.newLine(beamPoints[i-1].tag,beamPoints[i].tag)
        beamLines.append(l)

#Points and lines tendon
tendonSet=sets.defSet('tendonSet')
tendonPoints= tendonSet.getPoints
tendonLines= tendonSet.getLines
for i in range(1,3):
    p= points.newPntFromPos3d(geom.Pos3d(0,i*span/3,-deltaTendon))
    tendonPoints.append(p)
tendonLines.append(lines.newLine(beamPoints[0].tag,tendonPoints[0].tag))
tendonLines.append(lines.newLine(tendonPoints[0].tag,tendonPoints[1].tag))
tendonLines.append(lines.newLine(tendonPoints[1].tag,beamPoints[3].tag))


#BEAM
#Geometric section
from materials.sections import section_properties as sectpr
geomSectBeam=sectpr.RectangularSection(name='geomSectBeam',b=wBeam,h=hBeam)

# Material definition
concrete=tm.MaterialData(name='concrete',E=Ec,nu=nuc,rho=densc)
beamMat=tm.BeamMaterialData(name= 'beamMat', section=geomSectBeam, material=concrete)
beamMat.setupElasticShear3DSection(prep)

#Meshing
for l in beamSet.getLines:
    l.nDiv=nDivLines
beam_mesh=fem.LinSetToMesh(linSet=beamSet,matSect=beamMat,elemSize=None,vDirLAxZ=xc.Vector([1,0,0]),elemType='ElasticBeam3d',dimElemSpace=3,coordTransfType='linear')
beam_mesh.generateMesh(prep)

#Boundary conditions
modelSpace.fixNode000_FFF(0)
endnode=beamSet.getNodes.getNearestNode(geom.Pos3d(0,span,0))
modelSpace.fixNode000_FFF(endnode.tag)

#TENDON
#Material
prestressingSteel= tm.defCableMaterial(preprocessor=prep, name="prestressingSteel",E=Ep,prestress=fpi,rho=0.0)

#Meshing
for l in tendonSet.getLines:
    l.nDiv=nDivLines
corCooTr=modelSpace.newLinearCrdTransf(trfName='corCooTr',xzVector=xc.Vector([1,0,0]))
tendon_mesh=fem.LinSetToMesh(linSet=tendonSet,matSect=prestressingSteel,elemSize=None,vDirLAxZ=xc.Vector([1,0,0]),elemType='Truss',dimElemSpace=3,coordTransfType=None)
tendon_mesh.generateMesh(prep)
for e in tendonSet.getElements:
    e.area=Aps

# Connection between tendon and beam
gluedDOFs= [0,1,2,3,4,5]
for n1,n2 in zip(beamLines[0].getNodes(),tendonLines[0].getNodes()):
    modelSpace.constraints.newEqualDOF(n1.tag,n2.tag,xc.ID(gluedDOFs))
for n1,n2 in zip(beamLines[1].getNodes(),tendonLines[1].getNodes()):
    modelSpace.constraints.newEqualDOF(n1.tag,n2.tag,xc.ID(gluedDOFs))
for n1,n2 in zip(beamLines[2].getNodes(),tendonLines[2].getNodes()):
    modelSpace.constraints.newEqualDOF(n1.tag,n2.tag,xc.ID(gluedDOFs))

'''
#Plot
from postprocess.xcVtk.FE_model import vtk_FE_graphic
defDisplay= vtk_FE_graphic.RecordDefDisplayEF()
defDisplay.displayMesh(xcSets=[beamSet,tendonSet],fName= None,caption='Mesh',nodeSize=0.0010,scaleConstr=0.30)
'''
# Loads definition
cargas= prep.getLoadHandler
lPatterns= cargas.getLoadPatterns
#Load modulation.
ts= lPatterns.newTimeSeries("constant_ts","ts")
lPatterns.currentTimeSeries= "ts"
#Load case definition
lp0= lPatterns.newLoadPattern("default","0")
lPatterns.currentLoadPattern='0'
# #We add the load case to domain.
lPatterns.addToDomain('0') # THE ERROR IS HERE (WE ADD TO DOMAIN TOO EARLY)
print 'AAAAAAAAAAAAAAAAAAAA'

strain=(fps-fpi)/Ep
for e in tendonSet.getElements:
    eLoad= lp0.newElementalLoad("truss_temp_load")
    eLoad.elementTags= xc.ID([e.tag])
    eLoad.eps1= strain
    eLoad.eps2= strain
print 'BBBBBBBBBBBBBBBBBBBBBBBBBB'
#We add the load case to domain.
#lPatterns.addToDomain('0') #THE SOLUTION IS HERE

analisis= predefined_solutions.simple_static_linear(FEcase)
analOk= analisis.analyze(1)
print 'CCCCCCCCCCCCCCCCCCCCCCCCCCCCC'

from postprocess.xcVtk.FE_model import quick_graphics as QGrph
from postprocess.xcVtk import vtk_graphic_base
lcs=QGrph.QuickGraphics(FEcase)
lcs.displayDispRot(itemToDisp='uZ',setToDisplay=beamSet+tendonSet,fConvUnits=1e3,unitDescription='beam [mm]. Phase 1: prestressing of tendon',viewDef= vtk_graphic_base.CameraParameters("XNeg",1),fileName=None,defFScale=2e2)

lcs.displayIntForcDiag(itemToDisp='N',setToDisplay=tendonSet,fConvUnits= 1,scaleFactor=1,unitDescription='[kN,m]',viewDef=vtk_graphic_base.CameraParameters("ZNeg",1),fileName=None,defFScale=1)



