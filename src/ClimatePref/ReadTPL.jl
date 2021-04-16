using JuliaDB
tpldict = Dict(:ID => String,	:Major_group => String,	:Family => String,	:Genus_hybrid_marker => String,	:Genus => String,	:Species_hybrid_marker => String,	:Species => String,	:Infraspecific_rank => String,	:Infraspecific_epithet => String,	:Authorship => String,	:Taxonomic_status_in_TPL => String,	:Nomenclatural_status_original => String,	:Confidence_level => String,	:Source => String,	:Source_id => String,	:IPNI_id => String,	:Publication => String,	:Collation => String,	:Page => String,	:Date => String)
tplnames = ["ID",	"Major_group",	"Family",	"Genus_hybrid_marker",	"Genus",	"Species_hybrid_marker",	"Species",	"Infraspecific_rank",	"Infraspecific_epithet",	"Authorship",	"Taxonomic_status_in_TPL",	"Nomenclatural_status_original",	"Confidence_level",	"Source",	"Source_id",	"IPNI_id",	"Publication",	"Collation",	"Page",	"Date"]
function ReadTPL(folder::String)
    tpl = loadtable(folder, header_exists = false, colnames = tplnames, colparsers = tpldict, chunks = 1)
    tpl = filter(t-> t.Taxonomic_status_in_TPL == "Accepted" && t.Infraspecific_rank == "", tpl)
    duplicates = ["kew-2623303", "kew-2777599"]
    tpl = filter(t -> t.ID ∉ duplicates, tpl)
    spp = collect(select(tpl, :Species))
    gen = collect(select(tpl, :Genus))
    Species_name = gen .* " " .* spp
    tpl = pushcol(tpl, :species, Species_name)
    tpl = pushcol(tpl, :SppID, 1:length(tpl))
    tpl = reindex(tpl, :species)
    tpl = popcol(tpl, :Species)
    return tpl
end
