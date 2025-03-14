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

# # Manipulating Fields of floats
#
# Fields in MEDCoupling have a unique mesh support with a fixed and well-defined
# dimension.
#
# _(This may seem trivial, but it's actually a major difference from the concept
# of fields in the `MED-file` library, which is much more permissive.)_
#
# Fields are useful for:
#
# - storing values of a physical quantity over a domain
# - using some functions which compute values related to a given mesh, such as the
#   volume of its cells. This includes functions like:
#   - `getValueOn()`
#   - `getValueOnMulti()`
#   - `integral()`
#   - `getMeasureField`
#   - `normL1()`
#   - `normL2()`
#   - `fillFromAnalytic()`
# - precisely specifying the information exchanged between different codes during
#   during a coupling of those codes.
#
# For your information, the implementation of `MEDCouplingFieldDouble` is
# relatively small because this class the vast majority of its processing to other
# underlying classes like `MEDCouplingMesh`, `DataArrayDouble`, and
# `MEDCouplingSpatialDiscretization`. The `MEDCouplingFieldDouble` class ensures
# consistency among all these elements.
#
# It's often possible, and sometimes even recommended, to directly manipulate the
# arrays (a `DataArrayDouble`) and/or the mesh of an instance of
# `MEDCouplingFieldDouble`.
#

# ## Objectives
#
# This exercise focuses on the relationship between the meshes and the values of a
# field.
#
# - Creating a field
# - Aggregating fields
# - Building a subset of a field
# - Renumbering entities of a field
# - Comparing two fields from different sources
# - Evaluating a field at a set of points
# - Exploding a field
#

#
# We import the Python module `medcoupling`.
#

import medcoupling as mc

# We create a `MEDCouplingUMesh` from a 3D Cartesian mesh. Each direction will contain 10 cells and 11 nodes.

xarr = mc.DataArrayDouble.New(11, 1)
xarr.iota(0.0)  # Generate s, s+1, s+2, ... given start value s
cmesh = mc.MEDCouplingCMesh.New()
cmesh.setCoords(xarr, xarr, xarr)
mesh = cmesh.buildUnstructured()

# **Note:** The method `MEDCouplingMesh.buildUnstructured()` is very useful for
# quickly constructing an unstructured mesh to test something.  To highlight the
# problem of multiple geometric types, we convert the cells with even identifiers
# to polyhedra.

mesh.convertToPolyTypes(mc.DataArrayInt.Range(0, mesh.getNumberOfCells(), 2))

# ## Creating a field
#
# Create a scalar field (a single component) on cells (i.e. a "P0" field, in
# the MEDCoupling vocabulary), from the following analytical
# function: `f: x -> ( x - 5) * ( x - 5) + (y - 5) * (y - 5) + (z - 5) * (z - 5)`, where `(x, y, z)`
# represents the coordinates of the centroid of a cell. Two possibilities:
#
# - Directly by calling `fillFromAnalytic()` on a mesh
#

f = mesh.fillFromAnalytic(
    mc.ON_CELLS, 1, "(x-5.)*(x-5.)+(y-5.)*(y-5.)+(z-5.)*(z-5.)"
)  # 1 means that the field should have one component
f.setName("MyField")

# - Or by first creating an uninitialized field and applying `fillFromAnalytic()`
#   on this instance of `MEDCouplingFieldDouble`
#

f2 = mc.MEDCouplingFieldDouble(mc.ON_CELLS, mc.ONE_TIME)
f2.setMesh(mesh)
f2.setName("MyField2")
f2.fillFromAnalytic(
    1, "(x-5.)*(x-5.)+(y-5.)*(y-5.)+(z-5.)*(z-5.)"
)  # 1 means that the field should have one component

# Compare the two fields: compare `f` and `f2` with a precision of 1e-12 on the coordinates and 1e-13 on the values.

print("Are f and f2 equal?", f.isEqualWithoutConsideringStr(f2, 1e-12, 1e-13))

# <div class="alert alert-block alert-success">
# <b>Note:</b> The "WithoutConsideringStr" in the name of the previous method indicates that the field names will not be compared. This suffix can be found in other MEDCoupling methods.</div>
#
# ## Building a subset of a field
#
# Retrieve in a variable `ids1` the list of cell identifiers for which the field value is in the range [0.0,5.0]. Use the method `DataArrayDouble.findIdsInRange()`. With this result, build the subset `fPart1` of field `f`.

# this DataArrayDouble, which is a direct reference (not a copy) of the field's values
da1 = f.getArray()
ids1 = da1.findIdsInRange(0.0, 5.0)
fPart1 = f.buildSubPart(ids1)
fPart1.writeVTK("ExoField_fPart1.vtu")

# <img src="FieldDouble1.png" style="width:500px;">
#
# Select the part `fPart2` of field `f` where all tuple values are in `[50.,+infinity)`.

ids2 = f.getArray().findIdsInRange(50.0, 1.0e300)
fPart2 = f.buildSubPart(ids2)

# This kind of technique makes it easier to extract parts of a field related to a group of cells, for example.
#
# ## Renumbering entities of a field
#
# The generated `fPart1` part is valid from a MEDCoupling perspective. But it's not valid from a MED file perspective. Renumbering is necessary if you intend to store this field in a MED file to order cells by geometric type.
#
# The idea is to use the two methods `MEDCouplingUMesh.sortCellsInMEDFileFrmt()` and `DataArrayDouble.renumberInPlace()` to manually renumber a copy of `fPart1`:

fPart1Cpy = fPart1.deepCopy()
o2n = fPart1Cpy.getMesh().sortCellsInMEDFileFrmt()
fPart1Cpy.getArray().renumberInPlace(o2n)

# `fPart1Cpy` is now normalized to be stored in a MED file (which we'll see later).
#
# Check that `fPart1Cpy` and `fPart1` are the same apart from a permutation (`MEDCouplingFieldDouble.substractInPlaceDM()`)

fPart1Cpy.substractInPlaceDM(fPart1, 12, 1e-12)
fPart1Cpy.getArray().abs()
print("Equal field ? %s" % (fPart1Cpy.getArray().accumulate()[0] < 1e-12))

# <div class="alert alert-block alert
#
# -success">
# <b>Note:</b> The renumbering performed here is actually a very particular case of interpolation. Indeed, the assumption is made that the supports of fPart1 and fPart1Cpy are equal up to a permutation of cells and/or nodes.</div>
#
# ## Aggregating fields
#
# Aggregate `fPart1` and `fPart2` (use `MEDCouplingFieldDouble.MergeFields()`). And put the result of the aggregation in `fPart12`.

fPart12 = mc.MEDCouplingFieldDouble.MergeFields([fPart1, fPart2])
fPart12.writeVTK("ExoField_fPart12.vtu")

# <div class="alert alert-block alert-success">
# <b>Note:</b> The method MEDCouplingFieldDouble.MergeFields() should really be named MEDCouplingFieldDouble.AggregateFields() ...</div>
#
# <img src="FieldDouble2.png" style="width:500px;">
#
# ## Evaluating a field at given points in space
#
# Evaluate the value of field `fPart12` calculated previously at the cell centroids of its mesh (variable `bary`) and put the result in `arr1`. Use the methods `MEDCouplingFieldDouble.getValueOnMulti()` and `MEDCouplingMesh.computeCellCenterOfMass()` for this.
#
# Similarly, then evaluate the field `f` directly using the same list of points as before (`bary`) and put the result in `arr2`.
#
# Then check that `arr1` and `arr2` are indeed equal:

bary = fPart12.getMesh().computeCellCenterOfMass()
arr1 = fPart12.getValueOnMulti(bary)
arr2 = f.getValueOnMulti(bary)
delta = arr1 - arr2
delta.abs()
print("Is field evaluation matching?", (delta.accumulate()[0] < 1e-12))

# <div class="alert alert-block alert-success">
# <b>Note:</b> In this context, and for a cell-centered (P0) field for example, "evaluating" at a point means returning the value of the cell containing the given point. For node-centered (P1) fields, cells must be of simple types (triangles, tetrahedra) and linear interpolation is then used.</div>
#
# <div class="alert alert-block alert-success">
# <b>Note:</b> This technique can be used to quickly assess the quality of an interpolation.</div>
#
# ## Operations on fields
#
# Compute the integral of field `fPart12` over the mesh, and find it in another way using the method `DataArrayDouble.accumulate()` on the value array of this field. Recall that, given the simplified mesh in play, all cells have a unit volume.

integ1 = fPart12.integral(0, True)
integ1_bis = fPart12.getArray().accumulate()[0]
print("First integral matching ?", (abs(integ1 - integ1_bis) < 1e-8))

# Then apply a homothety of factor 1.2 centered at [0.,0.,0.] on the support of `fPart12` (i.e., its mesh). What is the new value of the integral?

fPart12.getMesh().scale([0.0, 0.0, 0.0], 1.2)
integ2 = fPart12.integral(0, True)
print("Second integral matching ?", (abs(integ2 - integ1_bis * 1.2 * 1.2 * 1.2) < 1e-8))

# ## Exploding a field - Displacement vectors
#
# We will now create a new mesh representing the exploded mesh of the initial mesh.
#
# Starting from the `mesh` mesh, create a vector field on cells `fVec` having 3 components representing the displacement vector between the point [5.,5.,5.] and the centroid of each cell of the mesh. Use the method `MEDCouplingMesh.fillFromAnalytic()`:

fVec = mesh.fillFromAnalytic(mc.ON_CELLS, 3, "(x-5.)*IVec+(y-5.)*JVec+(z-5.)*KVec")

# <div class="alert alert-block alert-success">
# <b>Note:</b> The special identifiers IVec, JVec, and KVec represent the unit vectors of the frame.</div>
#
# Then create a reduction of `fVec` (named `fVecPart1`) on the cells `ids1` previously obtained:

fVecPart1 = fVec.buildSubPart(ids1)
fVecPart1.setName("fVecPart1")

# Build the scalar field `fPart1Exploded` having the same values as `fPart1` but based on an exploded mesh compared to that of `fPart1.getMesh()`. To explode `fPart1.getMesh()` use the vector displacement field `fVecPart1` to apply to each cell the associated translation.

# +
cells = fPart1.getMesh().getNumberOfCells() * [None]

for icell, vec in enumerate(fVecPart1.getArray()):
    m = fPart1.getMesh()[[icell]]
    m.zipCoords()  # Not mandatory but saves memory
    m.translate(vec)
    cells[icell] = m
    pass

meshFVecPart1Exploded = mc.MEDCouplingUMesh.MergeUMeshes(cells)
fPart1.setMesh(meshFVecPart1Exploded)
fPart1.writeVTK("ExoField_fPart1_explo.vtu")
# -

# And here's what should be obtained:
#
# <img src="FieldDouble3.png" style="width:500px;">
