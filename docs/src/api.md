```@meta
CurrentModule = EcoSISTEM
```

# API

```@index
```

**`EcoSISTEMRasterDataSourcesExt` is in this list, and no other extension is.** Every public name
whose implementation moved into an extension keeps its docstring on a method-less stub in the parent,
which is why the rest need no entry - but `Base.read`'s dataset methods **cannot** take a stub
(the parent would have to define `read(::Type, ...)`, pirating `Base.read` for every type), so their
docstrings have nowhere to live but the extension. Naming the module is what keeps them in the
manual.

```@autodocs
Modules = filter(!isnothing,
                 [EcoSISTEM, EcoSISTEM.ClimatePref, EcoSISTEM.Units,
                  Base.get_extension(EcoSISTEM,
                                     :EcoSISTEMRasterDataSourcesExt)])
Private = true
```

```@index
```
