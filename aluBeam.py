# -*- coding: utf-8 -*-
from __future__ import division 
from __future__ import print_function


import math
import xc_base
import geom
import xc

from model import predefined_spaces
from model.geometry import grid_model as gm
from model.mesh import finit_el_model as fem
from model.sets import sets_mng as sets
from materials import typical_materials as tm
from materials.astm_aisc import ASTM_materials
from actions import loads
from actions import load_cases as lcases
from postprocess import limit_state_data as lsd
from materials.sections.structural_shapes import aisc_metric_shapes as lb
from actions import combinations as combs
from materials.astm_aisc import AISC_limit_state_checking as aisc
from connections.steel_connections import cbfem_bolt_weld as cbfem

# Default configuration of environment variables.
from postprocess.config import default_config
from postprocess import output_styles as outSty
from postprocess import output_handler as outHndl


# ***Data***
#steelPlate=ASTM_materials.A36   #steel shear plates
steelBolt=ASTM_materials.A325   
# aluminium mechanical porperties
E_alu=70e9 #Young modulus (Pa)
fty_alu=250e6 # yield strength (Pa)
ftu_alu=290e6 # ultimate strength (Pa)
nu_alu=0.33  #poisson coefficient
dens_alu=2700 # density (kg/m3)
# Geometry
height=150e-3    #beam height
length=0.5  #length of beam modelized
flng_width=80e-3 #flange width
flng_thck=3e-3   #flange thickness

web_thck=2e-3   #web thickness
web_int_dist=8e-3 #internal distance between webs
tab_int_dist=4.8e-3 #internal distance between tabs at the top of the profile
bot_alveolus_int_height=44e-3 #internal height of the bottom alveolus between webs
mid_alveolus_int_height=43.2e-3 #internal height of the middle alveolus between webs
top_alveolus_int_height=44e-3 #internal height of the top alveolus between webs
stiff_thck=web_thck #thickness of the bottom and middle stiffeners between webs

bolt_diam=8e-3  # diameter of bolt (only used to generate an instance of Bolt class)
hole_diam=10.3e-3 #diameter of the holes
hole_xdist=40e-3 # distance between holes in transverse direction
hole_ydist=30e-3 # distance between holes in longitudinal direction
hole_ydist_edge=25e-3 # distance from the first hole to the edge
nhole_yrows= 3 # number of rows of holes (longitudinal direction)

# *** derived parameters
hole_rad=hole_diam/2
tab_thck=(web_int_dist+2*web_thck)/2 #thickness of the top tabs
topstiff_thck=tab_thck #thickness of the top stiffener between webs
z_bottom_flange=0
z_top_flange=height-flng_thck
z_bot_stiff=flng_thck/2+bot_alveolus_int_height+stiff_thck/2
z_mid_stiff=z_bot_stiff+mid_alveolus_int_height+stiff_thck
z_top_stiff=z_mid_stiff+top_alveolus_int_height+topstiff_thck

x_web=web_int_dist/2+web_thck/2
x_flange=flng_width/2

x_holes=[-hole_xdist/2,hole_xdist/2]
y_holes=[hole_ydist_edge+i*hole_ydist_edge for i in range(nhole_yrows)]
esize=20e-3
angRadii=math.radians(45)

F=2e3
# **********End data and parameters ***************

xLst=[-x_flange,-x_web,0,x_web,x_flange] ; xLst.sort()  
yLst=[0,length] ; yLst.sort() 
#yLst=[0,-yb,yb,yCp] ; yLst.sort() 
zLst=[z_bottom_flange,z_bot_stiff,z_mid_stiff,z_top_stiff,z_top_flange] ; zLst.sort() 

FEcase= xc.FEProblem()
FEcase.title= 'Beam domo alusphere'
preprocessor=FEcase.getPreprocessor
prep=preprocessor   #short name
nodes= prep.getNodeHandler
elements= prep.getElementHandler
elements.dimElem= 3
# Problem type
# dimension of the space: nodes by three coordinates (x,y,z) and 
# six DOF for each node (Ux,Uy,Uz,thetaX,thetaY,thetaZ)
modelSpace= predefined_spaces.StructuralMechanics3D(nodes) #Defines the
surfaces= modelSpace.getSurfaceHandler()
points=modelSpace.getPointHandler()

sty=outSty.OutputStyle() 
out=outHndl.OutputHandler(modelSpace,sty)

grid=gm.GridModel(prep,xLst,yLst,zLst)
grid.generatePoints()

botFlange=grid.genSurfOneXYZRegion(((-x_flange,0,z_bottom_flange),(x_flange,length,z_bottom_flange)),'botFlange')

topFlange=grid.genSurfMultiXYZRegion([[(-x_flange,0,z_top_flange),(-x_web,length,z_top_flange)],
                                      [(x_web,0,z_top_flange),(x_flange,length,z_top_flange)]],'topFlange')

web=grid.genSurfMultiXYZRegion([[(-x_web,0,z_bottom_flange),(-x_web,length,z_top_stiff)],
                                [(x_web,0,z_bottom_flange),(x_web,length,z_top_stiff)]],'web')

tabs=grid.genSurfMultiXYZRegion([[(-x_web,0,z_top_stiff),(-x_web,length,z_top_flange)],
                                [(x_web,0,z_top_stiff),(x_web,length,z_top_flange)]],'tabs')

lowStiff=grid.genSurfMultiXYZRegion([[(-x_web,0,z_bot_stiff),(x_web,length,z_bot_stiff)],
                                     [(-x_web,0,z_mid_stiff),(x_web,length,z_mid_stiff)]],'lowStiff')

topStiff=grid.genSurfOneXYZRegion(((-x_web,0,z_top_stiff),(x_web,length,z_top_stiff)),'topStiff')

lstSurfSets=[botFlange,topFlange,web,tabs,lowStiff,topStiff]
for st in lstSurfSets:
    for s in st.surfaces:
        s.setElemSize(esize,True)



bolt=cbfem.Bolt(bolt_diam,steelBolt)

def gen_bolt_radii_pnts(x,ylist,z):
    # list of lists: [point_of_bolt,[8 points of radii]]
    lst_bolts_radii=list()
    for y in ylist:
        lst_radii=list()
        pBolt=points.newPntFromPos3d(geom.Pos3d(x,y,z))
        pHole=[points.newPntFromPos3d(geom.Pos3d(x+hole_rad*math.cos(i*angRadii),y+hole_rad*math.sin(i*angRadii),z)) for i in range(8)]
        pTags=[p.tag for p in pHole]
        pFace=surfaces.newPolygonalFacePts(pTags)
        lst_bolts_radii.append([pBolt,pHole,pFace])
    return lst_bolts_radii


negTopFlangeSet=grid.getSetSurfOneXYZRegion(((-x_flange,0,z_top_flange),(-x_web,length,z_top_flange)),'negTopFlangeSet')

negTopFlange=negTopFlangeSet.surfaces[0]
negTopFlangeBoltPnt= gen_bolt_radii_pnts(-hole_xdist/2,y_holes,z_top_flange) #list of lists: [point_of_bolt,[8 points of hole]]
for blt in negTopFlangeBoltPnt:
    h=blt[2]
    h.setNDiv(1)
    negTopFlange.addHole(h)
#meshSet.surfaces.append(negTopFlange)


posTopFlangeSet=grid.getSetSurfOneXYZRegion(((x_web,0,z_top_flange),(x_flange,length,z_top_flange)),'posTopFlangeSet')
posTopFlange=posTopFlangeSet.surfaces[0]
posTopFlangeBoltPnt= gen_bolt_radii_pnts(hole_xdist/2,y_holes,z_top_flange) #list of lists: [point_of_bolt,[8 points of hole]]
for blt in posTopFlangeBoltPnt:
    h=blt[2]
    h.setNDiv(1)
    posTopFlange.addHole(h)
#meshSet.surfaces.append(posTopFlange)


negBotFlangeSet=grid.getSetSurfOneXYZRegion(((-x_flange,0,z_bottom_flange),(-x_web,length,z_bottom_flange)),'negBotFlangeSet')
negBotFlange=negBotFlangeSet.surfaces[0]
negBotFlangeBoltPnt= gen_bolt_radii_pnts(-hole_xdist/2,y_holes,z_bottom_flange) #list of lists: [point_of_bolt,[8 points of hole]] 
for blt in negBotFlangeBoltPnt:
    h=blt[2]
    h.setNDiv(1)
    negBotFlange.addHole(h)
#meshSet.surfaces.append(negBotFlange)

posBotFlangeSet=grid.getSetSurfOneXYZRegion(((x_web,0,z_bottom_flange),(x_flange,length,z_bottom_flange)),'posBotFlangeSet')
posBotFlange=posBotFlangeSet.surfaces[0]
posBotFlangeBoltPnt= gen_bolt_radii_pnts(hole_xdist/2,y_holes,z_bottom_flange) #list of lists: [point_of_bolt,[8 points of hole]]
for blt in posBotFlangeBoltPnt:
    h=blt[2]
    h.setNDiv(1)
    posBotFlange.addHole(h)
#meshSet.surfaces.append(posBotFlange)

centBotFlangeSet=grid.getSetSurfOneXYZRegion(((-x_web,0,z_bottom_flange),(x_web,length,z_bottom_flange)),'centBotFlangeSet')

### Define material
#flange_mat= tm.defElasticMembranePlateSection(preprocessor, "flange_mat",E=E_alu,nu=nu_alu,rho= dens_alu,h= flng_thck)
# Materials for linear analysis vonmises verification
aluminium=tm.defElasticIsotropic3d(preprocessor=prep, name='aluminium', E=E_alu, nu=nu_alu, rho=dens_alu)
flange_mat=tm.DeckMaterialData(name='flange_mat',thickness= flng_thck,material=aluminium)
flange_mat.setupElasticSection(preprocessor=prep)
### Define template element
seedElemHandler= preprocessor.getElementHandler.seedElemHandler
seedElemHandler.defaultMaterial= flange_mat.name
elem= seedElemHandler.newElement("ShellMITC4",xc.ID([0,0,0,0]))

xcTotalSet= modelSpace.getTotalSet()
# Rename surfaces to avoid error in getTag()
for i, s in enumerate(xcTotalSet.surfaces):
    s.name= 'f'+str(i)

    
    
meshSet= modelSpace.defSet('meshSet')

for s in [negTopFlange,posTopFlange,negBotFlange,posBotFlange]:
    meshSet.surfaces.append(s)

for st in [web,tabs,lowStiff,topStiff,centBotFlangeSet]:
    for s in st.surfaces:
        meshSet.surfaces.append(s)

meshSet.useGmsh= True
meshSet.fillDownwards()
meshSet.genMesh(xc.meshDir.I)

negTopFlangeSet.fillDownwards()


for lstBoltHoles in [negTopFlangeBoltPnt,posTopFlangeBoltPnt,negBotFlangeBoltPnt,posBotFlangeBoltPnt]:
    for blt in lstBoltHoles:
        pntBolt=blt[0]
        pntBolt.genMesh(xc.meshDir.I) 
        bolt.generateRadii(tPlate=flng_thck,pntBolt=pntBolt,pntHoleLst=blt[1])
        nBolt=pntBolt.getNode()
        modelSpace.fixNode000_FFF(nBolt.tag)
    
out.displayFEMesh([topFlange,botFlange,web,tabs,lowStiff,topStiff])


# Assign thickness to elements
for e in topFlange.elements:
    e.getPhysicalProperties.getVectorMaterials[0].h=flng_thck
    e.setProp('yieldStress',fty_alu)

for e in botFlange.elements:
    e.getPhysicalProperties.getVectorMaterials[0].h=flng_thck
    e.setProp('yieldStress',fty_alu)

for e in web.elements:
    e.getPhysicalProperties.getVectorMaterials[0].h=web_thck
    e.setProp('yieldStress',fty_alu)

for e in tabs.elements:
    e.getPhysicalProperties.getVectorMaterials[0].h=tab_thck
    e.setProp('yieldStress',fty_alu)

for e in lowStiff.elements:
    e.getPhysicalProperties.getVectorMaterials[0].h=stiff_thck
    e.setProp('yieldStress',fty_alu)
          
for e in topStiff.elements:
    e.getPhysicalProperties.getVectorMaterials[0].h=topstiff_thck
    e.setProp('yieldStress',fty_alu)
         

# Analysis
loadNodes=sets.get_set_nodes_plane_XZ(setName='loadNodes',setBusq=xcTotalSet,yCoord=length, tol=0.0001)   
nnodes=loadNodes.nodes.size
Fnode=F/nnodes
lstLoadNod=[n for n in loadNodes.nodes]
nodalLoad
