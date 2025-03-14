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

# + [markdown] editable=true slideshow={"slide_type": ""}
# # Reading, Writing a MED file using the advanced MEDLoader API
#
# The advanced MEDLoader API is represented by the `MEDFile*` classes of the MEDLoader library.
#
# > - At the highest level, for the entire file: `MEDFileData`,
# > - For all meshes in the file: `MEDFileMeshes`,
# > - For each mesh: `MEDFileMeshMultiTS`, `MEDFileMesh`, `MEDFileUMesh`, `MEDFileCMesh`,
# > - For all fields in the file: `MEDFileFields`, `MEDFileFieldGlobs`,
# > - And finally, for each field: `MEDFileField1TS`, `MEDFileFieldMultiTS`
#
# ## Objective
#
# Write a mesh and a field from scratch, reread them, and compare the results.
#
# Topics covered: using the advanced MEDLoader API
#
# > - Writing a file
# > - Reading a file
#
# ## Implementation Start
#
# This exercise, like all others, relies on the Python scripting language. We load the Python module `MEDLoader`.
#
# For your information, the complete `MEDCoupling` module is included in `MEDLoader`. There's no need to import it if `MEDLoader` has been loaded.
# -

import medcoupling as mc

# + [markdown] editable=true slideshow={"slide_type": ""}
# ## Mesh Reading, Writing
#
# First, let's create the same mesh `targetMesh` as for the simple API.

# + editable=false slideshow={"slide_type": ""}
# fmt: off
targetCoords = [
    -0.3, -0.3, 0.2, -0.3, 0.7, -0.3, -0.3, 0.2, 0.2,
     0.2,  0.7, 0.2, -0.3, 0.7,  0.2,  0.7, 0.7, 0.7,
]
# fmt: on
qua1 = [0, 3, 4, 1]
tri1 = [1, 4, 2]
tri2 = [4, 5, 2]
qua2 = [6, 7, 4, 3]
qua3 = [7, 8, 5, 4]
targetConnTest = qua1 + tri1 + tri2 + qua2 + qua3
# -

targetMesh = mc.MEDCouplingUMesh("MyMesh", 2)
targetMesh.allocateCells(5)
targetMesh.insertNextCell(mc.NORM_TRI3, 3, tri1)
targetMesh.insertNextCell(mc.NORM_TRI3, 3, tri2)
targetMesh.insertNextCell(mc.NORM_QUAD4, 4, qua1)
targetMesh.insertNextCell(mc.NORM_QUAD4, 4, qua2)
targetMesh.insertNextCell(mc.NORM_QUAD4, 4, qua3)
myCoords = mc.DataArrayDouble(targetCoords, 9, 2)
myCoords.setInfoOnComponents(["X [km]", "YY [mm]"])
targetMesh.setCoords(myCoords)

# <div class="alert alert-block alert-success">
# <b>Note:</b> The mesh `targetMesh` is ordered by geometric type.</div>
#
# Next, we construct `targetMesh1`, representing the subcomponents (faces) of mesh `targetMesh`, and extract only the cells (here, surfaces) [3,4,7,8]. For more details on descending connectivity, refer to the **Descending Connectivity** section of the second exercise. This set could, for example, represent an area of interest for a calculation:

targetMeshConsti, _, _, _, _ = targetMesh.buildDescendingConnectivity()
targetMesh1 = targetMeshConsti[[3, 4, 7, 8]]
targetMesh1.setName(targetMesh.getName())

# <div class="alert alert-block alert-success">
# <b>Note:</b> In Python, the underscore `_` indicates that we expect a return value but won't use it (we don't bind it).</div>
#
# <div class="alert alert-block alert-success">
# <b>Note:</b> `targetMesh1` will be saved as part of the same global mesh in the MED file. Therefore, it must have the same name. This illustrates how a mesh in the MED file sense can mix dimensions.</div>
#
# Then, we write both meshes to the file "TargetMesh2.med".

meshMEDFile = mc.MEDFileUMesh()
meshMEDFile.setMeshAtLevel(0, targetMesh)
meshMEDFile.setMeshAtLevel(-1, targetMesh1)
meshMEDFile.write("TargetMesh2.med", 2)  # 2 stands for 'write from scratch'

# ## Reading, Writing Cell Groups
#
# Let's create two cell groups on the 2D mesh, i.e., at the relative level 0 (here, relative level 0 corresponds to 2D, level -1 corresponds to 1D, etc.). The first group `grp0_Lev0` contains cells [0,1,3], the second `grp1_Lev0` contains cells [1,2,3,4]:

grp0_0 = mc.DataArrayInt([0, 1, 3])
grp0_0.setName("grp0_Lev0")
grp1_0 = mc.DataArrayInt([1, 2, 3, 4])
grp1_0.setName("grp1_Lev0")
meshMEDFile.setGroupsAtLevel(0, [grp0_0, grp1_0])

# <div class="alert alert-block alert-success">
# <b>Note:</b> Naming the arrays is crucial; the name will be used for the group.</div>
#
# Let's create three level -1 groups, i.e., face groups. The first one called `grp0_LevM1` contains cells [0,1], the second one called `grp1_LevM1` contains cells [0,1,2], and the third one `grp2_LevM1` contains cells [1,2,3]:

grp0_M1 = mc.DataArrayInt([0, 1])
grp0_M1.setName("grp0_LevM1")
grp1_M1 = mc.DataArrayInt([0, 1, 2])
grp1_M1.setName("grp1_LevM1")
grp2_M1 = mc.DataArrayInt([1, 2, 3])
grp2_M1.setName("grp2_LevM1")
meshMEDFile.setGroupsAtLevel(-1, [grp0_M1, grp1_M1, grp2_M1])

# Write it all:

meshMEDFile.write("TargetMesh2.med", 2)  # 2 stands for 'write from scratch'

# Then, we can reread the MED file:

meshMEDFileRead = mc.MEDFileMesh.New(
    "TargetMesh2.med"
)  # a new is needed because it returns a MEDFileUMesh (MEDFileMesh is abstract)
meshRead0 = meshMEDFileRead.getMeshAtLevel(0)
meshRead1 = meshMEDFileRead.getMeshAtLevel(-1)
print(
    "Is level 0 in the file equal to 'targetMesh'?",
    meshRead0.isEqual(targetMesh, 1e-12),
)
print(
    "Is level 0 in the file equal to 'targetMesh1'?",
    meshRead1.isEqual(targetMesh1, 1e-12),
)

# Display the available levels for the group `grp0_Lev0`:

print(meshMEDFileRead.getGrpNonEmptyLevels("grp0_Lev0"))

# And finally, retrieve the cell identifiers contained in the group `grp0_Lev0`:

grp0_0_read = meshMEDFileRead.getGroupArr(0, "grp0_Lev0")
print(
    "Is group 'grp0_Lev0' equal to what is read in the file?",
    grp0_0_read.isEqual(grp0_0),
)

# ## Reading/Writing Fields with the Advanced API
#
# Let's create a simple vector field, on cells (P0), with a single timestep, called `f`.

f = mc.MEDCouplingFieldDouble(mc.ON_CELLS, mc.ONE_TIME)
f.setTime(5.6, 7, 8)
f.setArray(targetMesh.computeCellCenterOfMass())
f.setMesh(targetMesh)
f.setName("AFieldName")

# Store `f` in an object `MEDFileField1TS` (a field with a single timestep) to prepare for MED writing.

fMEDFile = mc.MEDFileField1TS()
fMEDFile.setFieldNoProfileSBT(f)  # No profile desired on the field, Sort By Type

# Add the field to the file "TargetMesh2.med".

fMEDFile.write(
    "TargetMesh2.med", 0
)  # 0 is paramount to indicate that we *append* (and no overwrite) to the MED file

# <div class="alert alert-block alert-success">
# <b>Note:</b> Note the use of 0 to indicate that we want to add to the existing file.</div>
#
# Read the field:

fMEDFileRead = mc.MEDFileField1TS("TargetMesh2.med", f.getName(), 7, 8)
fRead1 = fMEDFileRead.getFieldOnMeshAtLevel(
    mc.ON_CELLS, 0, meshMEDFileRead
)  # Quickest way, not re-reading mesh in the file.
fRead2 = fMEDFileRead.getFieldAtLevel(
    mc.ON_CELLS, 0
)  # Like above, but this time the mesh is read!
print(
    "Does the field remain OK with the quick method?", fRead1.isEqual(f, 1e-12, 1e-12)
)
print("Does the field remain OK with the slow method?", fRead2.isEqual(f, 1e-12, 1e-12))

# ## Reading/Writing a Field on a Profile
#
# Now, let's see an advanced concept of MED files, the ability to write a field on only a part of the mesh. The commonly used technique is to put specific values (e.g., +infinity, i.e., 1e+300) on areas where the field doesn't make sense, thus helping to detect any bugs during the calculation.
#
# The operation mode with profiles remains uncommon.
#
# Let's construct a reduction to cells [1,2,3] of `f` and call it `fPart`:

pfl = mc.DataArrayInt([1, 2, 3])
pfl.setName("My1stPfl")
fPart = f.buildSubPart(pfl)
fPart.setName("fPart")

# Store it in the `MEDFileField1TS` structure and invoke `setFieldProfile()`.

fMEDFile2 = mc.MEDFileField1TS()
fMEDFile2.setFieldProfile(
    fPart, meshMEDFileRead, 0, pfl
)  # 0 is the relative level (here 0 means 2D)
fMEDFile2.write(
    "TargetMesh2.med", 0
)  # 0 is paramount to indicate that we *append* (and no overwrite) to the MED file

# Read the `fPart` field from the file "TargetMesh2.med" and the corresponding cell identifiers.

fMEDFileRead2 = mc.MEDFileField1TS("TargetMesh2.med", fPart.getName(), 7, 8)
fPartRead, pflRead = fMEDFileRead2.getFieldWithProfile(mc.ON_CELLS, 0, meshMEDFileRead)
print(
    "Is the partial field correctly read?",
    fPartRead.isEqualWithoutConsideringStr(fPart.getArray(), 1e-12),
)
print(
    "Is the list of cell identifiers matching?",
    pflRead.isEqualWithoutConsideringStr(pfl),
)
