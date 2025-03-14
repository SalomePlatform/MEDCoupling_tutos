# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.16.7
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# # MEDCoupling / MEDLoader - Exemple complet 1 - Agitateur
#
# Here, we start with a file called agitateur.med with the following content:
#
# <img src="agitateur.jpg" style="width:600px;">
#
# This is the result of a small two-phase calculation: the green magnetic stirrer (identified only by a cell field and having no proper mesh) rotates from one time step to another within a liquid phase. Two drops of liquid fall during this time towards the air/water interface (in gray).
#
# The goal of the exercise is to calculate the torque applied to this stirrer, which is the mechanical piece driving the lower part of the fluid.
#
# ## Objective
#
# The objective is to provide a complete example of non-trivial post-processing from a MED file.
#
# ## Start of the implementation
#
# To start the exercise, import the entire python module `MEDLoader` (which includes `MEDCoupling`). Also import `numpy`.

import medcoupling as mc
import numpy as np

# ## **Extracting meshes and fields with the advanced API**
#
# With the advanced API, read the entire file `agitateur.med` and display all the time steps of the 1st field.

# +
filename: str = "agitateur.med"
data = mc.MEDFileData(filename)
f0 = data.getFields()[0]
print(type(f0))
print(f"f0 name: {f0.getName()}")

ts = f0.getTimeSteps()
print(f"ts: {ts}")
# -

# Retrieve the mesh of the stirrer (in green) at time step (2, -1) (cf. ts). The position of the stirrer is defined by a field on the global mesh of the system and has no proper mesh. Therefore, use the cell field `DISTANCE_INTERFACE_ELEM_BODY_ELEM_DOM` and select only the part of the field with a value within `[0.,1.]`. Put the corresponding cell identifiers in `ids`:

fMts = data.getFields()["DISTANCE_INTERFACE_ELEM_BODY_ELEM_DOM"]
f1ts = fMts[(2, -1)]
fMc = f1ts.getFieldAtLevel(mc.ON_CELLS, 0)
arr = fMc.getArray()
arr.getMinMaxPerComponent()  # just to see the field variation range per component
ids = arr.findIdsInRange(0.0, 1.0)
f2Mc = fMc[ids]

# Using the field named `PRESSION_ELEM_DOM`, we find the 3D pressure field applied by the stirrer.

pressMts = data.getFields()["PRESSION_ELEM_DOM"]
press1ts = pressMts[(2, -1)]
pressMc = press1ts.getFieldAtLevel(mc.ON_CELLS, 0)
pressOnAgitatorMc = pressMc[ids]

# Remove unnecessary nodes from `pressOnAgitatorMc.getMesh()`:

pressOnAgitatorMc.getMesh().zipCoords()

# ## **Transition from a 3D cell field to a 3D surface field**
#
# Deduct the 3D pressure field on the surface of the stirrer. To do this, pass through the descending mesh `MEDCouplingUMesh.buildDescendingConnectivity()`.

agitateurMesh3DMc = pressOnAgitatorMc.getMesh()
(
    m3DSurf,
    desc,
    descI,
    revDesc,
    revDescI,
) = agitateurMesh3DMc.buildDescendingConnectivity()
nbOf3DCellSharing = revDescI.deltaShiftIndex()
ids2 = nbOf3DCellSharing.findIdsEqual(
    1
)  # Cells with only one neighbor are on the boundary, i.e. on the skin
agitateurSkinMc = m3DSurf[ids2]
offsetsOfTupleIdsInField = revDescI[ids2]
tupleIdsInField = revDesc[offsetsOfTupleIdsInField]
pressOnSkinAgitateurMc = pressOnAgitatorMc[tupleIdsInField]
pressOnSkinAgitateurMc.setMesh(agitateurSkinMc)

# ## **Manipulating fields**
#
# Calculate the vector field of force on the stirrer's surface by multiplying the pressure for each cell by the surface and then by the normal vector. The pressure is in bar, so convert it beforehand to pascal (Pa).

pressSkin = pressOnSkinAgitateurMc.getArray()
pressSkin *= 1e5  # conversion from bar to Pa
areaSkin = agitateurSkinMc.getMeasureField(True).getArray()
forceSkin = pressSkin * areaSkin
normalSkin = agitateurSkinMc.buildOrthogonalField().getArray()
forceVectSkin = forceSkin * normalSkin

# Now, let's calculate the first moment at the center of mass of the stirrer:
#
# To make this first calculation of the torque exerted on the stirrer, let's calculate the position of the stirrer's center of mass. Calculate the polyhedron representing the envelope of the 3D mesh of the stirrer `agitateurMesh3DMc` (use `MEDCouplingUMesh.buildSpreadZonesWithPoly()`).

singlePolyhedron = agitateurMesh3DMc.buildSpreadZonesWithPoly()
singlePolyhedron.orientCorrectlyPolyhedrons()
centerOfMass = singlePolyhedron.computeCellCenterOfMass()

# <div class="alert alert-block alert-success">
# <b>Note:</b> Calling MEDCouplingUMesh.orientCorrectlyPolyhedrons() is not obligatory but recommended because if by chance the polyhedron is incorrectly oriented, its barycenter will be incorrect!</div>
#
# Calculate for each cell of the stirrer's surface the moment with respect to the center of mass `centerOfMass` of the stirrer. To do this, calculate `posSkin`, the `DataArrayDouble` giving for each cell of the stirrer's surface the vector `centerOfMass` -> `G`, with `G` the barycenter of the current cell.

barySkin = agitateurSkinMc.computeCellCenterOfMass()
posSkin = barySkin - centerOfMass

# Now apply the classical formula for moment calculation: calculate the cross product by cell of `posSkin` with `forceVectSkin` (method `DataArrayDouble.CrossProduct()`).

torquePerCellOnSkin = mc.DataArrayDouble.CrossProduct(posSkin, forceVectSkin)

# Sum `torqueOnSkin` using the method `DataArrayDouble.accumulate()`.

zeTorque = torquePerCellOnSkin.accumulate()
print("torque = %r N.m" % zeTorque[2])

# Let's check the previously calculated torque by dividing the power by the angular velocity. The linear velocity is stored in the field "VITESSE_ELEM_DOM".
#
# Calculate the power for each cell of the stirrer's surface and sum it.

speedMts = data.getFields()["VITESSE_ELEM_DOM"]
speed1ts = speedMts[(2, -1)]
speedMc = speed1ts.getFieldAtLevel(mc.ON_CELLS, 0)
speedOnSkin = speedMc.getArray()[tupleIdsInField]
powerSkin = mc.DataArrayDouble.Dot(forceVectSkin, speedOnSkin)
power = powerSkin.accumulate()[0]
print("power = %r W" % (power))

# Calculate the angular velocity. To do this, calculate the sum of `x^2`, `y^2`, and `xz` from `posSkin` and build (with NumPy) the 2x2 inertia matrix `inertiaSkin=[[x2,xy], [xy,z2]]`.
#
# Retrieve the eigenvector associated with the maximum eigenvalue with `linalg.eig(inertiaSkin)`.

x2 = posSkin[:, 0] * posSkin[:, 0]
x2 = x2.accumulate()[0]
y2 = posSkin[:, 1] * posSkin[:, 1]
y2 = y2.accumulate()[0]
xy = posSkin[:, 0] * posSkin[:, 1]
xy = xy.accumulate()[0]
inertiaSkin = np.matrix([[x2, xy], [xy, y2]])
inertiaSkinValues, inertiaSkinVects = np.linalg.eig(inertiaSkin)
pos = max(enumerate(inertiaSkinValues), key=lambda x: x[1])[0]
vect0 = inertiaSkinVects[pos].tolist()[0]
print(vect0)

# With the previous calculation, we can deduce that the stirrer has rotated by 1.1183827931 radians (cf. complete solution for details - we put the previous steps in a function that we apply to several time steps).
#
# Calculate and compare the torque on the stirrer.

omega = 1.1183827931 / (ts[-1][2] - ts[0][2])
print(
    "At timestep (%d,%d) (physical time=%r s) the torque is: %r N.m, power/omega=%r N.m "
    % (ts[2][0], ts[2][1], ts[2][2], zeTorque[2], power / omega)
)
