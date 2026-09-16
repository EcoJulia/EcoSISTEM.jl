# SPDX-License-Identifier: LGPL-3.0-or-later

# Taking a built layer apart into plain data and putting it back together.
#
# A layer is split into a **skeleton** - its concrete type, its dims, its cell size and its change,
# with every array replaced by a descriptor - and the **arrays** themselves, in the order the rebuild
# consumes them. The skeleton holds no function and no array values, so it is small and can be
# serialised; the arrays can be sent in place by anything that moves plain `Array`s. This is how the
# MPI extension hands one rank's build to the others.
#
# Only what `materialise` builds from a spec is handled: `ContinuousLayer` and `CategoricalLayer`
# over a `DimArray`, with no change or an absolute `SeriesLayerChange`, and a `LayerCollection` of
# those; and a `Raster`, which is what a study area's `active` is. Anything else is refused by name,
# so a new layer or change type fails here until it is taught.

# The largest array a descriptor may describe: MPI counts elements in a C `int`.
const _PAYLOADLIMIT = Int(typemax(Int32))

# Take `layer` apart. Returns the skeleton, the arrays in rebuild order, and one descriptor per array
# saying what storage a receiver must allocate for it.
function _layerpayload(layer::AbstractLayer)
    arrays = Array[]
    descriptors = NamedTuple[]
    skeleton = _packlayer!(arrays, descriptors, layer, "the layer")
    return (skeleton = skeleton, arrays = arrays, descriptors = descriptors)
end

# Put a layer back together from `_layerpayload`'s skeleton and arrays. The arrays are consumed in
# order, and every one of them must be.
function _rebuildlayer(skeleton::NamedTuple, arrays::AbstractVector)
    queue = collect(Array, arrays)
    layer = _unpacklayer(skeleton.type, skeleton, queue)
    isempty(queue) ||
        error("a layer payload carried $(length(queue)) more array(s) than its skeleton describes.")
    return layer
end

# Uninitialised storage for one described array, for a receiver to fill.
function _allocatepayload(descriptor::NamedTuple)
    return Array{descriptor.eltype}(undef, descriptor.size)
end

# Refuse an array MPI cannot send in one piece, before anything is sent or allocated.
function _checkdescriptor(descriptor::NamedTuple, label::AbstractString)
    n = prod(descriptor.size)
    n <= _PAYLOADLIMIT ||
        error("$label has $n elements, and an array is sent between MPI ranks in one piece, which " *
              "allows at most $_PAYLOADLIMIT. Use a coarser grid or a shorter series.")
    return descriptor
end

# One layer's skeleton. The two leaf types hold the same three fields, and `typeof(layer)` carries
# every type parameter, so rebuilding by it restores the role, axis, eltype and array type.
function _packlayer!(arrays, descriptors,
                     layer::Union{ContinuousLayer, CategoricalLayer}, label)
    return (type = typeof(layer),
            matrix = _packdimarray!(arrays, descriptors, layer.matrix,
                                    "$label's `matrix`"),
            size = layer.size,
            change = _packchange!(arrays, descriptors, layer.change,
                                  "$label's change"))
end

# A collection's members, in name order, with the names kept.
function _packlayer!(arrays, descriptors, coll::LayerCollection, label)
    nt = getfield(coll, :nt)
    return (type = typeof(coll), names = keys(nt),
            members = Tuple(_packlayer!(arrays, descriptors, member,
                                        "layer `$name`")
                            for (name, member) in pairs(nt)))
end

# Any other layer type is refused, naming it.
function _packlayer!(arrays, descriptors, layer::AbstractLayer, label)
    return error("$label is a $(nameof(typeof(layer))), which cannot yet be sent between MPI " *
                 "ranks.")
end

# No change: nothing but its type.
function _packchange!(arrays, descriptors, change::NoLayerChange, label)
    return (type = NoLayerChange,)
end

# A read series: the times, origin and policies go in the skeleton; the slices and the baseline
# are arrays.
function _packchange!(arrays, descriptors,
                      change::SeriesLayerChange{AbsoluteChange}, label)
    return (type = typeof(change),
            slices = _packdata!(arrays, descriptors, change.slices,
                                "$label's `slices`"),
            times = change.times, origin = change.origin,
            atend = change.atend, calendar = change.calendar,
            baseline = _packdimarray!(arrays, descriptors, change.baseline,
                                      "$label's `baseline`"))
end

# Any other change is refused, naming it - a declared change may hold a function.
function _packchange!(arrays, descriptors, change::AbstractLayerChange, label)
    return error("$label is a $(nameof(typeof(change))), which cannot yet be sent between MPI " *
                 "ranks.")
end

# A `DimArray`: dims, reference dims, name and metadata in the skeleton, the values as an array.
function _packdimarray!(arrays, descriptors, x::DimArray, label)
    return (type = typeof(x), dims = dims(x),
            refdims = DimensionalData.refdims(x),
            name = DimensionalData.name(x),
            metadata = DimensionalData.metadata(x),
            data = _packdata!(arrays, descriptors, parent(x), label))
end

# A `Raster` - a study area's `active` is one - is a `DimArray` with a missing value as well.
function _packdimarray!(arrays, descriptors, x::Rasters.Raster, label)
    return (type = typeof(x), dims = dims(x),
            refdims = DimensionalData.refdims(x),
            name = DimensionalData.name(x),
            metadata = DimensionalData.metadata(x),
            missingval = Rasters.missingval(x),
            data = _packdata!(arrays, descriptors, parent(x), label))
end

# Any other array type is refused, naming it.
function _packdimarray!(arrays, descriptors, x, label)
    return error("$label is a $(nameof(typeof(x))), which cannot yet be sent between MPI ranks.")
end

# The values themselves. A `BitArray` is sent as `Bool`s, since its packed chunks are not the array.
function _packdata!(arrays, descriptors, data::Array, label)
    descriptor = _checkdescriptor((eltype = eltype(data), size = size(data),
                                   bits = false), label)
    push!(arrays, data)
    push!(descriptors, descriptor)
    return descriptor
end

function _packdata!(arrays, descriptors, data::BitArray, label)
    descriptor = _checkdescriptor((eltype = Bool, size = size(data),
                                   bits = true), label)
    push!(arrays, Array{Bool}(data))
    push!(descriptors, descriptor)
    return descriptor
end

# Values held in anything but a plain array (a view, say) are refused, naming the type.
function _packdata!(arrays, descriptors, data, label)
    return error("$label is held in a $(nameof(typeof(data))), which cannot yet be sent between " *
                 "MPI ranks.")
end

# The inverses of `_packlayer!`, dispatched on the type the skeleton recorded.
function _unpacklayer(::Type{T}, skeleton,
                      queue) where {T <: Union{ContinuousLayer,
                                          CategoricalLayer}}
    matrix = _unpackdimarray(skeleton.matrix, queue)
    return T(matrix, skeleton.size,
             _unpackchange(skeleton.change.type, skeleton.change, queue))
end

function _unpacklayer(::Type{<:LayerCollection}, skeleton, queue)
    members = map(m -> _unpacklayer(m.type, m, queue), skeleton.members)
    return LayerCollection(NamedTuple{skeleton.names}(members))
end

# The inverses of `_packchange!`, likewise.
_unpackchange(::Type{NoLayerChange}, skeleton, queue) = NoLayerChange()

function _unpackchange(::Type{T}, skeleton,
                       queue) where {T <: SeriesLayerChange}
    slices = _unpackdata(skeleton.slices, queue)
    baseline = _unpackdimarray(skeleton.baseline, queue)
    return T(slices, skeleton.times, skeleton.origin, skeleton.atend,
             skeleton.calendar, baseline)
end

# Rebuilt through the constructor rather than `rebuild`, since a receiver has no template; the
# result must come out as the very type that was sent.
function _unpackdimarray(skeleton, queue)
    x = _newdimarray(skeleton.type, _unpackdata(skeleton.data, queue),
                     skeleton)
    typeof(x) === skeleton.type ||
        error("a `DimArray` was rebuilt as a $(typeof(x)), not the $(skeleton.type) that was sent.")
    return x
end

# The next array in the queue, checked against its descriptor and restored to bits if it was sent
# as `Bool`s.
function _unpackdata(descriptor, queue)
    isempty(queue) &&
        error("a layer payload ran out of arrays before its skeleton was complete.")
    data = popfirst!(queue)
    (eltype(data) === descriptor.eltype && size(data) == descriptor.size) ||
        error("a layer payload's array is $(eltype(data)) of size $(size(data)), where its " *
              "skeleton describes $(descriptor.eltype) of size $(descriptor.size).")
    return descriptor.bits ? BitArray(data) : data
end

# A `DimArray` or a `Raster` from its values and its skeleton.
function _newdimarray(::Type{<:DimArray}, data, skeleton)
    return DimArray(data, skeleton.dims, refdims = skeleton.refdims,
                    name = skeleton.name, metadata = skeleton.metadata)
end

function _newdimarray(::Type{<:Rasters.Raster}, data, skeleton)
    return Rasters.Raster(data, skeleton.dims, refdims = skeleton.refdims,
                          name = skeleton.name, metadata = skeleton.metadata,
                          missingval = skeleton.missingval)
end
