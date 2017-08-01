
function TraitFun(eco::Ecosystem, pos::Int64, spp::Int64)
  TraitFun(getenvtype(eco), eco, pos, spp)
end

function TraitFun(env::Type{ContinuousHab{Temp}}, eco::Ecosystem, pos::Int64, spp::Int64)
  T = gethabitat(eco, pos)
  T_opt, Var = getpref(env, eco, spp)
  return gettraitfun(eco)(T, T_opt, Var)
end

function TraitFun(env::Type{DiscreteHab{None}}, eco::Ecosystem, pos::Int64, spp::Int64)
  currentniche = gethabitat(eco, pos)
  preference = getpref(env, eco, spp)
  return gettraitfun(eco)(currentniche, preference)
end

function TraitFun(env::Type{ContinuousHab{None}}, eco::Ecosystem, pos::Int64, spp::Int64)
  return 1.0
end

function gettraitfun(eco::Ecosystem)
  return eco.relationship.traitfun
end

function getpref{E <: Envtype}(env::E, eco::Ecosystem, spp::Int64)
  getpref(env, eco, spp)
end

function getpref(env::Type{ContinuousHab{Temp}}, eco::Ecosystem, spp::Int64)
  traits = eco.spplist.traits
  return traits.mean[spp], traits.var[spp]
end

function getpref(env::Type{DiscreteHab{None}}, eco::Ecosystem, spp::Int64)
  return eco.spplist.traits.trait[spp]
end
