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

# + [markdown] editable=true slideshow={"slide_type": "slide"}
# # Reading, writing a MED file using the basic MEDLoader API
#
# The basic MEDLoader API is contained in the `MEDLoader` class. All methods of this class are static (they do not depend on a particular instance of the class), their names begin with a capital letter. All read/write operations are executed with each method call, and no internal state of the class is saved.

# + [markdown] editable=false slideshow={"slide_type": "subslide"}
# ## Objective
#
# Write a mesh and a field from scratch, reread them, and compare the results.
#
# Topics covered: using the basic `MEDLoader` API
#
# > - Writing a file
# > - Reading a file

# + [markdown] editable=false slideshow={"slide_type": "subslide"}
# ## Implementation Start
#
# This exercise, like all others, relies on the Python scripting language. We load the Python module `MEDLoader`.
#
# For your information, the complete `MEDCoupling` module is included in `MEDLoader`. There's no need to import it if `MEDLoader` has been loaded.

# + editable=false slideshow={"slide_type": "skip"}
import medcoupling as mc

# + [markdown] editable=false slideshow={"slide_type": "slide"}
# ## Mesh Reading, Writing
#
# First, let's create a mesh `targetMesh` composed of several geometric types.

# + editable=false slideshow={"slide_type": "subslide"}
# fmt: off
targetCoords = [
    -0.3, -0.3, 0.2, -0.3, 0.7, -0.3, -0.3, 0.2, 0.2,
     0.2,  0.7, 0.2, -0.3, 0.7,  0.2,  0.7, 0.7, 0.7,
]
# fmt: on
tri1 = [1, 4, 2]
tri2 = [4, 5, 2]
qua1 = [0, 3, 4, 1]
qua2 = [6, 7, 4, 3]
qua3 = [7, 8, 5, 4]

targetConn = tri1 + tri2 + qua1 + qua2 + qua3
print(targetConn)

# + editable=true slideshow={"slide_type": ""}
targetMesh = mc.MEDCouplingUMesh("MyMesh", 2)
targetMesh.allocateCells(5)
targetMesh.insertNextCell(mc.NORM_TRI3, 3, tri1)
targetMesh.insertNextCell(mc.NORM_TRI3, 3, tri2)
targetMesh.insertNextCell(mc.NORM_QUAD4, 4, qua1)
targetMesh.insertNextCell(mc.NORM_QUAD4, 4, qua2)
targetMesh.insertNextCell(mc.NORM_QUAD4, 4, qua3)
myCoords = mc.DataArrayDouble(targetCoords, 9, 2)
myCoords.setInfoOnComponents(["X [km]", "Y [mm]"])
targetMesh.setCoords(myCoords)

# + [markdown] editable=false slideshow={"slide_type": ""}
# <div class="alert alert-block alert-success">
# <b>Note:</b> The mesh targetMesh is ordered by geometric type.</div>
#
# The mesh can then be directly written...

# + editable=true slideshow={"slide_type": ""}
mc.WriteUMesh("TargetMesh.med", targetMesh, True)  # True means 'from scratch'

# + [markdown] editable=false slideshow={"slide_type": ""}
# ... and read.
# -

meshRead = mc.ReadUMeshFromFile("TargetMesh.med", targetMesh.getName(), 0)
print("Is the read mesh equal to 'targetMesh' ?", meshRead.isEqual(targetMesh, 1e-12))

# + [markdown] editable=false slideshow={"slide_type": ""}
# ## Read/Write a Field on a Time Step
#
# We now create a vector field `f` on cells (P0) with `targetMesh` as its support. This field corresponds, for example, to the physical time 5.6, marked by iteration 7 and sub-iteration 8. We take this opportunity to remind you that in `MEDCoupling` fields, physical time is given for information only; storage and most API functions are based on the last two integers.
# -

f = mc.MEDCouplingFieldDouble.New(mc.ON_CELLS, mc.ONE_TIME)
f.setTime(5.6, 7, 8)  # Declare the timestep associated with the field
f.setArray(targetMesh.computeCellCenterOfMass())
f.setMesh(targetMesh)
f.setName("AFieldName")
mc.WriteField("MyFirstField.med", f, True)

# + [markdown] editable=false slideshow={"slide_type": ""}
# Subsidiary question: What does the field created correspond to?
#
# <div class="alert alert-block alert-success">
# <b>Note:</b> The mesh and the field are written in one go to the file "MyFirstField.med".</div>
#
# We then read MyFirstField.med:
# -

f2 = mc.ReadFieldCell("MyFirstField.med", f.getMesh().getName(), 0, f.getName(), 7, 8)
print("Is the read field identical to 'f' ?", f2.isEqual(f, 1e-12, 1e-12))

# + [markdown] editable=false slideshow={"slide_type": ""}
# <div class="alert alert-block alert-success">
# <b>Note:</b> When reading the field, we must know its name, the name of its support mesh, and the desired timestep. Functions like MEDFileFields.getFieldsNames() or MEDFileMeshes.getMeshesNames() help with this.</div>
#
# <div class="alert alert-block alert-success">
# <b>Note:</b> The name ReadFieldCell() reminds us that the field must be read on the cells. Remember that according to the MED file standard, the same field can have some of its data stored on cells, but also simultaneously on nodes, Gauss points, etc... even if this kind of exotic mix is generally not recommended.</div>
#
# ## Read/Write a Field on Multiple Time Steps
#
# Here, unlike the previous case, we write multiple times to the same MED file. First, let's write the mesh.
# -

mc.WriteUMesh("MySecondField.med", f.getMesh(), True)

# + [markdown] editable=false slideshow={"slide_type": ""}
# Then, we write only the information related to the field (primarily its values array).
# -

mc.WriteFieldUsingAlreadyWrittenMesh("MySecondField.med", f)  # mesh is not re-written

# + [markdown] editable=false slideshow={"slide_type": ""}
# Next, we add a second timestep on the same mesh.
# -

f2 = f.clone(True)  # 'True' means that we need a deep copy
f2.getArray()[:] = 2.0
f2.setTime(7.8, 9, 10)
mc.WriteFieldUsingAlreadyWrittenMesh("MySecondField.med", f2)

# + [markdown] editable=false slideshow={"slide_type": ""}
# Now the file "MySecondField.med" contains the mesh and a field with two timesteps carried by this mesh.
#
# We can reread all of this with methods similar to what was seen previously:

# + editable=true slideshow={"slide_type": ""}
f3 = mc.ReadFieldCell("MySecondField.med", f.getMesh().getName(), 0, f.getName(), 7, 8)
print("Is the field read in the file equal to 'f' ?", f.isEqual(f3, 1e-12, 1e-12))
f4 = mc.ReadFieldCell("MySecondField.med", f.getMesh().getName(), 0, f.getName(), 9, 10)
print("Is the field read in the file equal to 'f2' ?", f2.isEqual(f4, 1e-12, 1e-12))
