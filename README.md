# MEDCoupling Tutorials

Welcome to the **MEDCoupling Tutorials** repository. This project provides
learning materials and practical examples to help you get started with
MEDCoupling and explore its capabilities.

These tutorials are designed to be used alongside the [MEDCoupling reference
manual](https://docs.salome-platform.org/latest/dev/MEDCoupling/developer/index.html),
which provides the full API documentation and detailed descriptions of
available methods.

## Repository Contents

This repository is organized as follows:

- **Course slides (source)**: `./presentations/src`
- **Course slides (PDF)**: `./presentations/install`
- **Jupyter notebook tutorials**: `./tutorials`
- **Visualization utilities**: `./src/medcoupling_tutos/`
  - Includes helper tools to convert MEDCoupling meshes and fields into
    `pyvista` objects for inline visualization in notebooks.

---

## Installation

You have two main options:

1. **Use MEDCoupling via Salome** and manually install the required Python
   dependencies (`numpy`, `jupyterlab`, `pyvista`, etc.).
2. **Use the simplified setup below** (recommended), which relies on `uv` for
   environment and dependency management.

---

## Quick Installation (Recommended)

### Linux

#### Requirements

Most distributions already include:

- `curl`
- `git`

#### Steps

```sh
# Install uv (fast Python package and environment manager)
curl -LsSf https://astral.sh/uv/install.sh | sh

# Clone the repository
git clone https://github.com/SalomePlatform/MEDCoupling_tutos.git
cd MEDCoupling_tutos

# Create environment and install dependencies
uv sync --python 3.12
```

---

### Windows

#### Requirements

- `git` (must be installed beforehand)

#### Steps

Install `uv`:

```powershell
powershell -ExecutionPolicy ByPass -c "irm https://astral.sh/uv/install.ps1 | iex"
```

Clone the repository and install dependencies:

```sh
git clone https://github.com/SalomePlatform/MEDCoupling_tutos.git
cd MEDCoupling_tutos

uv sync --python 3.12
```

---

## Launching the Tutorials

Start JupyterLab to explore the notebooks:

```sh
cd tutorials
uv run jupyter-lab
```

Then open the notebooks from your browser and follow the tutorials interactively.
You should use the [MEDCoupling
documentation](https://docs.salome-platform.org/latest/dev/MEDCoupling/developer/index.html)
alongside the tutorials.

---

## Contributing

Contributions are welcome! Feel free to open issues or submit pull requests to
improve tutorials, fix bugs, or extend functionality.
