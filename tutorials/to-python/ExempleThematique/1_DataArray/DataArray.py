# ## Manipulating `DataArray`
#
# ### Introduction
#
# `DataArray` (`DataArrayInt` and `DataArrayDouble`) are used in MEDCoupling to store values as contiguous arrays in memory. Values are grouped by tuples, and each tuple has the same number of components. They form the basis of many operations performed in MEDCoupling. Therefore, it is important to know how to manipulate them effectively.
#
# `DataArrayDouble` is often used for the direct manipulation of field values, as we will see later on.
#
# `DataArrayInt`, on the other hand, is used for the manipulation of cell and point identifiers.

#
# ### Summary
#
# The objective of this exercise is to become familiar with manipulating `DataArray`.
#
# The task is to create a mesh containing 4 adjacent squares.
#
# <img src="four_squares.svg" style="width:300px;">
#
# The first step consists in **creating a `DataArray`** (array), containing the coordinates of the nodes of a square, with a side length of $1$, centered at $(0, 0)$. The first component of the coordinate array is called `X` with the unit `m` (meter), and the 2nd component is called `Y` (same unit).
#
# Then, we construct the coordinates of the nodes of the four squares we want to obtain by **duplicating the coordinates** of the nodes of the original square, and then applying a **translation** to them.
#
# It is then necessary to **merge the duplicate nodes**, resulting from multiple duplications.
#
# Finally, an **unstructured mesh** is created from the `DataArray` of nodes which has been created.
#
# The concepts covered in this exercise are:
#
# - Creating an instance of `DataArrayDouble`
# - Displaying an instance of `DataArrayDouble` and invoking the `getValue()` method to convert it to a list
# - Using "slice" type notations such as `da[:,:]`
# - Learning about renumbering (old-2-new convention)
# - Invoking services such as `findCommonTuples()`
# - Creating a `UMesh` (unstructured mesh)

#
# ### Importing MEDCoupling
#
# It is necessary to import the Python module `medcoupling` to use its functionalities. We use the alias `mc` so that we don't have to rewrite `medcoupling` every time we call a function.
#
# All static methods of the module start with a capital letter. With these imports, the following are available:
#
# - all classes of MEDCoupling
# - all enumerations (for example, standard cell types: `mc.ON_CELLS`, `mc.ON_NODES`, `mc.ONE_TIME`...)
# - all static methods
#
# The native Python module `math` is also imported to perform trigonometric manipulations.

import medcoupling as mc
import math

# ### Creating an instance of `DataArrayDouble`, containing 4 tuples
#
# The goal here is to create an instance of `DataArrayDouble`, which will be used to store the coordinates of the points at the corners of a single square.
#
# There are different ways to construct a `DataArray`:

data_arr = mc.DataArrayDouble(4, 2)

# Here, the arguments `4` and `2` correspond to the **dimensions** of the `DataArray`.
#
# The first argument `4` indicates the number of **rows** of the array. This corresponds here to the number of corners of a square. For example, it would have been necessary to put `6` for a hexagon.
#
# The second argument `2` indicates the number of **columns**. This corresponds here to the number of coordinates for each point. Since we are working in a 2D space in this exercise, each point has two coordinates, relating to the $x$ and $y$ axes.

# Alternatively, we can first create a `DataArray` and specify its dimensions **later**.

# +
# allocate
data_arr = mc.DataArrayDouble()

# specify the dimension
data_arr.alloc(4, 2)

print(data_arr)
# -

# We can also construct a `DataArray` with **8 rows** and then rearrange it into **2 columns**.

# +
data_arr = mc.DataArrayDouble(8)

data_arr.rearrange(2)

print(data_arr)
# -

# **Note:** For the three approaches described above, it is important to note that the values of the array are **not initialized**. They contain meaningless values.

# Finally, we can construct a `DataArray` directly from a Python `list`. By default, the array has **only one component**. Therefore, it needs to be rearranged into **2 columns**.

data_arr = mc.DataArrayDouble(list(range(8)))
data_arr.rearrange(2)
print(data_arr)

# All these approaches are different ways to construct a `DataArrayDouble` containing 8 values, grouped into 4 tuples, each containing 2 components. Each tuple represents the 2D coordinates of a point.

# ### Displaying an array
#
# So far, the content of a `DataArray` has been displayed using the `print()` function, taking the array directly as an argument. This approach displays many meta-information about the object. It is also possible to display only the raw values of the array using the following syntax:

print(data_arr.getValues())

# ### Assigning values
#
# Assigning values to a `DataArray` is very similar to the syntax used for `ndarrays` from the `numpy` library. Notably, we can use **slicing**, using the colon character `:`.

data_arr[1:4, 0] = 3.0

# For example, we can directly assign values to the second column of the array using a list.

data_arr[:, 1] = [1.0, 2.0, 3.0, 4.0]

# We assign values to the array so that it contains the Cartesian coordinates of the nodes of a square, centered at `(0, 0)`, with a side length of `1`.

data_arr[:, 0] = [-0.5, 0.5, 0.5, -0.5]
data_arr[:, 1] = [-0.5, -0.5, 0.5, 0.5]
print(data_arr)

# **Note:** If necessary, it is possible to work in polar coordinates and perform conversions between polar and Cartesian coordinates.
#

# +
data_arr_pol = data_arr.fromPolarToCart()
print(data_arr_pol)

data_arr_cart = data_arr_pol.fromCartToPolar()
print(data_arr_cart)
# -

# We can assign **names to the different components (columns) of a `DataArray`**. This allows for clear indications on how to interpret the data in the array. Notably, we can specify the names of the coordinates (Cartesian, here), and the units (meter, here).

data_arr.setInfoOnComponents(mc.svec(["X [m]", "Y [m]"]))
print(data_arr)

# In the context of this exercise, this is not essential. However, other more advanced functions **require this information**.

# The four nodes of the square for which we constructed the coordinates belong to the circle centered at
#
#  `(0, 0)` with a radius of $\frac{\sqrt{2}}{2}$. We can verify that the norm (`magnitude()`) of each tuple is indeed equal to $\frac{\sqrt{2}}{2}$, with a tolerance of `1.e-12`, using the `isUniform()` method.

# +
print(data_arr.magnitude())

tolerance = 1.0e-12

print(
    "Is norm of each tuple equal to '0.5 * sqrt(2)' ?",
    data_arr.magnitude().isUniform(math.sqrt(2) / 2.0, tolerance),
)
# -

# ### Duplication and aggregation
#
# We now construct the list `translationToPerform`, which contains a list of vectors of size 2. This list of size 4 (4 squares) contains the different translations to be performed to obtain the coordinates of the nodes of the four squares we want to build.

translations = [
    [0.5, 0.5],
    [0.5, 1.5],
    [1.5, 0.5],
    [1.5, 1.5],
]

# We create the four arrays containing the translated coordinates of the nodes of the original square.

# list of 'DataArrayDouble'
squares = []
for translation in translations:
    # Adding a vector to a set of coordinates does a translation. translation
    # could have been a DataArrayDouble too.
    squares.append(data_arr + translation)
print(squares)

# ### Aggregation of arrays
#
# From this list of instances of `DataArrayDouble`, we build a single `DataArrayDouble`, the result of aggregating the instances one after the other. For this, we use the `Aggregate` method.

squares_da = mc.DataArrayDouble.Aggregate(squares)  # 'da' means 'data array', here
print(squares_da)

# We have thus constructed a `DataArrayDouble`, containing the coordinates of the nodes of the four squares we want to build. An important remark is that the order of aggregation of the `DataArray` is preserved.
#
# **Remarks:**
#
# - The same applies for the aggregation of meshes and fields, in order to facilitate access to and identification of data. This is, for example, an essential difference with the `MEDFile` model, as we will see later.
# - There is also the method `mc.DataArrayDouble.Meld(arr)`, which allows for the aggregation of two `DataArray` component by component, i.e., concatenating arrays column by column rather than row by row.

# ### Finding equal tuples
#
# In the `DataArray` constructed previously, some tuples are identical. It is necessary to detect the equal tuples to avoid creating duplicate nodes. To find equal tuples, we use `findCommonTuples()`. To compare two tuples of floating point numbers with each other, this function takes an absolute tolerance as an argument.

# **Note:**
# You can use `help(mc.DataArrayDouble.findCommonTuples)` to get information on how to use this function's interface.

# +
print("number of tuples: ", squares_da.getNumberOfTuples())

# help(mc.DataArrayDouble.findCommonTuples)

absolute_tolerance = 1.0e-12

common_tuples, common_tuples_indices = squares_da.findCommonTuples(absolute_tolerance)
print("common_tuples:\n", common_tuples.getValues())
print("common_tuples_indices:\n", common_tuples_indices.getValues())
# -

# The function `findCommonTuples` returns two `DataArrayInt`. The first contains groups of indices of common nodes, concatenated together. The second list indicates the starting position of each group in the first list.
#
# This return format in two `DataArrayInt` is quite common in MEDCoupling for performance reasons. It is called **indirect indexing**. This format is particularly used for functions manipulating unstructured meshes.
#
# **Attention: the last element of the second array points outside the first array**.
# This last index is always present and ensures that treatments such as slices presented just after are always valid, without needing to specify the last group.
#
# Here is a diagram illustrating how indirect indexing works, in a general way:
#
# <img src="IndirectIndex.jpg" style="width:700px;">
#
# The number of groups is therefore equal to the size of the second array, the one containing the indices, minus one.

c = common_tuples
ci = list(common_tuples_indices)
for i, (j, k) in enumerate(zip(ci[:-1], ci[1:])):
    print("group", i, ": ", c[int(j) : int(k)].getValues())

# For example, we can verify that tuples `1` and `8`, those in group `0`, are indeed equal:

grp_index = 0
grp = c[int(ci[grp_index]) : int(ci[grp_index + 1])]
for i in grp:
    print(squares_da[i].getValues())

# ### Calculating Group Sizes in Indirect Indexing
#
# In an indirect indexing, we can calculate the size of each group without directly manipulating the index array (i.e., the second array). For this purpose, the `DataArrayInt.deltaShiftIndex` function is utilized.

print(common_tuples_indices.deltaShiftIndex().getValues())

# ### Building an "old-to-new" Array
#
# Thanks to the previously performed duplicate node detection, it's possible to create a new table where duplicates have been eliminated.
#
# Mathematically, merging duplicate nodes into a single node is akin to performing a surjection from a starting space `X` (here with 16 nodes) to an arrival space `Y` with 9 nodes.
#
# <img src="SurjectionDataArray.png" style="width:250px;">
#
# In MEDCoupling, representing this surjection is done through a `DataArrayInt`, storing indices, following the MEDCoupling convention called **"old-to-new"**.
#
# **Note:** The "old-to-new" convention is preferred for bijections as well (e.g., a permutation).
#
# In this convention, each element at index $i$ in this array contains the new identifier of tuple $i$ from `X`, in the arrival set `Y`.
#
# For example, `[0, 1, 2, 2]`, in "old-to-new" convention, means that:
# - `X[0]` corresponds to `Y[0]`
# - `X[1]` corresponds to `Y[1]`
# - `X[2]` corresponds to `Y[2]`
# - `X[3]` corresponds to `Y[2]`
#
# We'll construct this table to extract a subset of the starting coordinates and keep only unique tuples (non-duplicates).
#
# The static method `DataArrayInt.ConvertIndexArrayToO2N()` allows us to switch from the storage mode of this surjection, represented by an indirect indexing, to the "old-to-new" format. The `O2N` suffix of the function stands for "OldToNew". This function also returns $Card(Y)$, which is the number of tuples after removing duplicates, here.

o2n, newNbOfTuples = mc.DataArrayInt.ConvertIndexArrayToO2N(
    len(squares_da), common_tuples, common_tuples_indices
)
print("old-to-new list: ", o2n.getValues())
print("new number of tuples: ", newNbOfTuples)

# In the final mesh, we indeed want to have only nine nodes (refer to the diagram of the four squares at the beginning of the exercise).
#
# Now we can build the array of unique tuples using `o2n` and `newNbOfTuples`. We use the `DataArrayDouble.renumberAndReduce()` function for this purpose.

squares_nodes_no_duplicate = squares_da.renumberAndReduce(o2n, newNbOfTuples)

# It's possible to translate all tuples in a `DataArray` at once by directly adding a list of the same dimension as the tuples in the array:

squares_nodes_no_duplicate += [1.0, 1.0]

# ### Building an Unstructured Mesh
#
# We can now create a 2D unstructured mesh using the constructed array to define the mesh nodes.

# +
dimension = 2

mesh = mc.MEDCouplingUMesh("FourSquares", dimension)
mesh.setCoords(squares_nodes_no_duplicate)

print("Mesh dimension is", mesh.getMeshDimension())
print("Spatial dimension is", mesh.getCoords().getNumberOfComponents())
# -

# Currently, the created mesh only contains a list of nodes. It doesn't have any cells yet.
#
# To create cells in this mesh, we need to create a group of node indices of a cell to add, then use the `insertNextCell` function. This group of node indices is called the "cell connectivity table".
#
# First, we need to allocate memory for the cells to create:

mesh.allocateCells(4)

# Finally, thanks to the previously constructed "old-to-new" table, we can create the connectivity table of each cell of the mesh, i.e., the four squares.

print(o2n.getValues())
for i in range(4):
    # iter through o2n by chunks of 4 nodes at a time
    square_node_indices = o2n[4 * i : 4 * (i + 1)]
    mesh.insertNextCell(mc.NORM_POLYGON, square_node_indices.getValues())

# We check that the mesh doesn't contain any anomalies.

mesh.checkConsistencyLight()

# It's always a good idea to call this method after constructing a mesh. It ensures there are no gross errors, particularly in the connectivity table of the mesh.
#
# To visually verify that the mesh is correct, we export it to a file in `.vtu` format, which can be viewed in `Paraview` or the `ParaVIS` module of the `Salome` platform.

mesh.writeVTK(mesh.getName() + ".vtu")

# **Note:** Here, we export the mesh in `.vtu` format and not `.med` because MEDCoupling does not include, by default, the `MED-file` library necessary for generating `.med` files.
