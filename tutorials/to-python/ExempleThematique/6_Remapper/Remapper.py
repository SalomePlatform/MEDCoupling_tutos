# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.16.7
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# # MEDCouplingRemapper: interpolation of fields
#
# Here, we will perform an interpolation between two meshes `srcMesh` and `trgMesh`. To emphasize some subtleties of interpolation, we take a particular case where `srcMesh` is a refined mesh of `trgMesh` (with some cells cut more finely).
#
# To start the exercise, import the `medcoupling` module and the `MEDCouplingRemapper` class from the module.

import medcoupling as mc
from medcoupling import MEDCouplingRemapper

# ## Create the target mesh
#
# Construct the unstructured mesh `trgMesh` from a 2D Cartesian mesh 10x10 starting at point `[0.,0.]` and having a step of 1 in both X and Y:

arr = mc.DataArrayDouble(11)
arr.iota(0)
trgMesh = mc.MEDCouplingCMesh()
trgMesh.setCoords(arr, arr)
trgMesh = trgMesh.buildUnstructured()

# ## Create the source mesh
#
# Create a mesh `srcMesh` from a 2D Cartesian mesh of 20x20 cells also starting at point `[0.,0.]` and having a step of 0.5 in both X and Y:

arr = mc.DataArrayDouble(21)
arr.iota(0)
arr *= 0.5
srcMesh = mc.MEDCouplingCMesh()
srcMesh.setCoords(arr, arr)
srcMesh = srcMesh.buildUnstructured()

# To make the exercise more interesting, triangulate the first 20 cells of `srcMesh` using `MEDCouplingUMesh.simplexize()` (2D simplices are triangles). Set the result to `srcMesh`.

tmp = srcMesh[:20]  # Extract a sub-part of srcMesh
tmp.simplexize(0)
srcMesh = mc.MEDCouplingUMesh.MergeUMeshes([tmp, srcMesh[20:]])

# Interpolate with MEDCouplingRemapper
#
# We recall that to project a field from one mesh to another, we must first prepare the interpolation matrix containing the projection ratios.
#
# Calculate the first part of the interpolation matrix from `srcMesh` (discretized at cells - P0) to `trgMesh` (also discretized at cells). To do this, invoke `MEDCouplingRemapper.prepare()` on an instance (`remap`) of `MEDCouplingRemapper`.

remap = MEDCouplingRemapper()
remap.prepare(srcMesh, trgMesh, "P0P0")

# Check that the matrix calculated by the method is correct in our trivial case. To do this, retrieve in `myMatrix` the internal matrix returned by `MEDCouplingRemapper.getCrudeMatrix()`. This gives for each cell of `trgMesh` the identifiers of cells of `srcMesh` with which it intersects, and the corresponding intersection area.
#
# Make sure that for each cell of `trgMesh`, the sum of the areas always equals 1.

myMatrix = remap.getCrudeMatrix()
print(myMatrix)
sumByRows = mc.DataArrayDouble(len(myMatrix))
for i, wIt in enumerate(sumByRows):
    su = 0.0
    for it in myMatrix[i]:
        su += myMatrix[i][it]
    wIt[0] = su
print("Is interpolation well prepared?", sumByRows.isUniform(1.0, 1e-12))

# Construct a field on cells "srcField" built from the following analytical formula: `7-sqrt((x-5.)*(x-5.)+(y-5.)*(y-5.))`:

srcField = mc.MEDCouplingFieldDouble(mc.ON_CELLS, mc.ONE_TIME)
srcField.setMesh(srcMesh)
srcField.fillFromAnalytic(1, "7-sqrt((x-5.)*(x-5.)+(y-5.)*(y-5.))")
srcField.getArray().setInfoOnComponent(0, "powercell [W]")

# Here is what this field looks like:
#
# <img src="Remapper1.png" style="width:500px;">
#
# Apply interpolation with `MEDCouplingRemapper.transferField()`:

# +
# remap.transferField(srcField, 1e300)
# -

# <div class="alert alert-block alert-success">
# <b>Note:</b> 1e300 is a default value. This value will be systematically assigned to any cell of trgField not intersecting any cell of srcMesh. Typically, users set an enormous value to spot what is often a bug. However, other users, from the perspective of parallel interpolation for example, set 0.</div>
#
# <div class="alert alert-block alert-success">
# <b>Note:</b> An exception is thrown because srcField has no defined nature. We will see later in the impact of this attribute on the final result.</div>
#
# Set the nature of `srcField` to `IntensiveMaximum`. This means that the field should be interpreted as intensive (such as temperature).

srcField.setNature(mc.IntensiveMaximum)
trgFieldCV = remap.transferField(srcField, 1e300)

# Check that with the `IntensiveMaximum` nature, the integral of the field is preserved. However, the sum over the cells (accumulation) is not preserved!

# +
integSource = srcField.integral(True)[0]
integTarget = trgFieldCV.integral(True)[0]
print("IntensiveMaximum -- integrals: %lf == %lf" % (integSource, integTarget))

accSource = srcField.getArray().accumulate()[0]
accTarget = trgFieldCV.getArray().accumulate()[0]
print("IntensiveMaximum -- sums: %lf != %lf" % (accSource, accTarget))
# -

# Now set the nature of `srcField` to `ExtensiveConservation`. The field should be interpreted as extensive (such as power or volume).

srcField.setNature(mc.ExtensiveConservation)
trgFieldI = remap.transferField(srcField, 1e300)

# Check that with the `ExtensiveConservation` nature, the integral of the field is not preserved. However, the sum over the cells is preserved.

# +
integSource = srcField.integral(True)[0]
integTarget = trgFieldI.integral(True)[0]
print("ExtensiveConservation -- integrals: %lf != %lf" % (integSource, integTarget))

accSource = srcField.getArray().accumulate()[0]
accTarget = trgFieldI.getArray().accumulate()[0]
print("ExtensiveConservation -- sums: %lf == %lf" % (accSource, accTarget))
# -

# Visualize the fields with ParaView, or by writing them to a file.
