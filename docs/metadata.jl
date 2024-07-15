# SPDX-License-Identifier: LGPL-3.0-or-later

using Pkg

# Update EcoSISTEM folder packages 
Pkg.activate(".")
Pkg.update()

# Update examples folder packages
if isdir("examples")
    if isfile("examples/Project.toml")
        Pkg.activate("examples")
        Pkg.update()
        Pkg.rm("EcoSISTEM")
        Pkg.develop("EcoSISTEM")
    end
end

# Update docs folder packages
Pkg.activate("docs")
Pkg.update()
Pkg.rm("EcoSISTEM")
Pkg.develop("EcoSISTEM")

# Reformat files in package
using JuliaFormatter
using EcoSISTEM
format(EcoSISTEM)

module Metadata

# Fix metadata
using Dates
using Git
using TOML
using JSON
using DataStructures
using HTTP
using YAML

function read_project()
    project_d = TOML.parsefile("Project.toml")
    project = OrderedDict{String, Any}()
    for key in [
        "name",
        "uuid",
        "author",
        "version",
        "license",
        "deps",
        "weakdeps",
        "extensions",
        "compat",
        "authors",
        "extras",
        "targets"
    ]
        if haskey(project_d, key)
            val = project_d[key]
            if val isa AbstractDict
                d = OrderedDict{String, Any}()
                for k2 in sort(collect(keys(val)))
                    d[k2] = val[k2]
                end
                project[key] = d
            else
                project[key] = val
            end
            delete!(project_d, key)
        end
    end
    for key in keys(project_d)
        project[key] = project_d[key]
    end

    return project
end

function get_person_from_orcid(orcid::String)
    url = "https://pub.orcid.org/v3.0/$orcid"
    headers = ["Accept" => "application/json"]
    response = HTTP.get(url, headers)

    if response.status == 200
        data = JSON.parse(String(response.body))
        name = data["person"]["name"]
        given_names = name["given-names"]["value"]
        family_name = name["family-name"]["value"]
        d = Dict("orcid" => orcid, "pid" => data["orcid-identifier"]["uri"],
                 "givenName" => given_names, "familyName" => family_name,
                 "full_name" => "$given_names $family_name")
        emails = data["person"]["emails"]["email"]
        if length(emails) ≥ 1
            d["email"] = emails[1]["email"]
            d["name_with_email"] = d["full_name"] * " <" * d["email"] * ">"
        else
            d["name_with_email"] = d["full_name"]
        end
        return d
    else
        return nothing
    end
end

function get_organisation_from_ror(ror::String)
    url = "https://api.ror.org/organizations/$ror"

    response = HTTP.get(url)

    if response.status == 200
        data = JSON.parse(String(response.body))
        d = Dict("name" => data["name"], "ror" => ror, "pid" => data["id"])
        return d
    else
        return nothing
    end
end

function get_first_release_date()
    project = read_project()
    package = project["name"]
    url = "https://raw.githubusercontent.com/JuliaRegistries/General/master/$(package[1])/$package/Versions.toml"
    headers = ["Accept" => "application/toml"]
    response = HTTP.get(url, headers)

    if response.status == 200
        data = TOML.parse(String(response.body))
        version = minimum(VersionNumber.(keys(data)))
        date = readchomp(`$(Git.git()) log -1 --format=%ad --date=format:%Y-%m-%d refs/tags/v$version`)
        return date
    end

    return nothing
end

function increase_patch()
    project = read_project()
    version = project["version"]
    v = VersionNumber(version)
    new_version = VersionNumber(v.major, v.minor, v.patch + 1)
    @info "Bumping patch version from $version to $new_version"
    project["version"] = string(new_version)
    open("Project.toml", "w") do io
        return TOML.print(io, project)
    end
    return crosswalk()
end

function increase_minor()
    project = read_project()
    version = project["version"]
    v = VersionNumber(version)
    new_version = VersionNumber(v.major, v.minor + 1, 0)
    @info "Bumping minor version from $version to $new_version"
    project["version"] = string(new_version)
    open("Project.toml", "w") do io
        return TOML.print(io, project)
    end
    return crosswalk()
end

function increase_major()
    v = VersionNumber(version)
    new_version = VersionNumber(v.major + 1, 0, 0)
    @info "Bumping major version from $version to $new_version"
    project["version"] = string(new_version)
    open("Project.toml", "w") do io
        return TOML.print(io, project)
    end
    return crosswalk()
end

function get_os_from_workflows()
    files = filter(isfile, readdir(".github/workflows", join = true))
    oses = Set{String}()
    for file in files
        jobs = YAML.load_file(file)["jobs"]
        for job in keys(jobs)
            os = jobs[job]["runs-on"]
            if occursin(r"\${{.*}}", os)
                k2s = split(replace(os, r"\${{ *([^ ]*) *}}" => s"\1"), ".")
                os = jobs[job]["strategy"]
                for k2 in k2s
                    os = os[k2]
                end
                if os isa String
                    push!(oses, os)
                else
                    for val in os
                        push!(oses, val)
                    end
                end
            else
                push!(oses, os)
            end
        end
    end

    platforms = Set{String}()
    for os in oses
        push!(platforms,
              replace(os, "ubuntu" => "Linux", "windows" => "Windows",
                      r"-.*" => ""))
    end

    return sort(collect(platforms))
end

function crosswalk()
    now = string(today())
    init = readchomp(`$(Git.git()) log --max-parents=0 --format=%ad --date=short -n 1`)
    tags = readlines(`$(Git.git()) tag -l --sort="version:refname"`)
    tag = maximum(VersionNumber.(tags))
    tag_date = readchomp(`$(Git.git()) log -1 --format=%ad --date=format:%Y-%m-%d refs/tags/v$tag`)
    branch = readchomp(`$(Git.git()) branch --show-current`)
    remotes = split(readchomp(`$(Git.git()) remote`), '\n')
    urls = String[]
    for remote in remotes
        push!(urls,
              replace(readchomp(`$(Git.git()) remote get-url $remote`),
                      r"\.git" => ""))
    end

    repos = replace.(urls, r"^.*/([^/]+)$" => s"\1")

    codemeta = isfile("codemeta.json") ?
               JSON.parsefile("codemeta.json", dicttype = OrderedDict) :
               OrderedDict{String, Any}()

    codemeta["@context"] = "https://w3id.org/codemeta/3.0"
    codemeta["type"] = "SoftwareSourceCode"
    codemeta["applicationCategory"] = "ecology"
    codemeta["programmingLanguage"] = "julia"
    codemeta["developmentStatus"] = "active"

    repo_index = "origin" ∈ remotes ?
                 repo_index = findfirst(==("origin"), remotes) : 1

    if haskey(codemeta, "codeRepository")
        cm_url = codemeta["codeRepository"]
        if cm_url ∉ urls
            @error "codemeta has wrong repo URL – $cm_url not in $urls – using $(remotes[repo_index])"
        else
            repo_index = findfirst(==(cm_url), urls)
        end
    end
    codemeta["codeRepository"] = urls[repo_index]

    if haskey(codemeta, "name")
        cm_name = codemeta["name"]
        if cm_name ≠ repos[repo_index]
            @error "codemeta has wrong repo repo name – $cm_name not $(repos[repo_index]) – fixing"
        end
    end
    codemeta["name"] = repos[repo_index]

    codemeta["issueTracker"] = urls[repo_index] * "/issues"

    readme = urls[repo_index] * "/blob/" * branch * "/README.md"
    cm_readme = get!(codemeta, "readme", readme)
    cm_readme == readme ||
        @info "README set to $cm_readme, not $readme"

    if haskey(codemeta, "buildInstructions")
        codemeta["buildInstructions"] == readme ||
            @info "README set to $(codemeta["buildInstructions"]), not $readme"
    else
        @warn "No build instructions"
    end

    cm_created = get(codemeta, "dateCreated", "")
    if cm_created ≠ init
        @warn "Fixing creation date to first git commit: $init"
        codemeta["dateCreated"] = init
    end

    platforms = get_os_from_workflows()
    cm_platforms = sort(string.(get(codemeta, "operatingSystem", String[])))
    if length(platforms) ≠ length(cm_platforms) ||
       any(platforms .≠ cm_platforms)
        @error "codemeta platforms do not match workflows ($cm_platforms ≠ $platforms), fixing"
        codemeta["operatingSystem"] = platforms
    end

    years = string(year(Date(init)))

    project = read_project()
    proj_version = VersionNumber(project["version"])

    cm_version = VersionNumber(get!(codemeta, "version", string(proj_version)))

    if proj_version == tag
        @debug "Still on latest release version: $tag"
        codemeta["dateModified"] = tag_date
        if cm_version ≠ tag
            @warn "Correcting codemeta tag version ($cm_version) to release tag ($tag)"
            cm_version = tag
            tag_year = string(year(Date(tag_date)))
            if tag_year ≠ years
                years = years * "-" * tag_year
            end
        else
            tag_year = string(year(Date(codemeta["dateModified"])))
            if tag_year ≠ years
                years = years * "-" * tag_year
            end
        end
    elseif proj_version > tag
        @info "Preparing for new release"
        codemeta["dateModified"] = now
        if cm_version ≠ proj_version
            @info "Updating codemeta tag version ($cm_version) to " *
                  "new release ($proj_version)"
            cm_version = proj_version
            this_year = string(year(Date(now)))
            if this_year ≠ years
                years = years * "-" * this_year
            end
        else
            tag_year = string(year(Date(codemeta["dateModified"])))
            if tag_year ≠ years
                years = years * "-" * tag_year
            end
        end
    else # Project version is lower than latest release!
        @error "Project.toml version is behind release ($proj_version < $tag), fixing."
        proj_version = tag
        cm_version = tag
        codemeta["dateModified"] = tag_date
        tag_year = string(year(Date(codemeta["dateModified"])))
        if tag_year ≠ years
            years = years * "-" * tag_year
        end
    end

    first_release_date = get_first_release_date()
    if !isnothing(first_release_date)
        if haskey(codemeta, "datePublished")
            codemeta["datePublished"] == first_release_date ||
                @warn "codemeta.json publication date inconsistent with Julia's General registry, fixing ($(codemeta["datePublished"]) ≠ $first_release_date)"
        end
        codemeta["datePublished"] = first_release_date
    end
    project["version"] = string(proj_version)
    codemeta["version"] = "v$cm_version"

    codemeta["downloadUrl"] = urls[repo_index] * "/archive/refs/tags/" *
                              codemeta["version"] * ".tar.gz"

    authors = String[]
    author_data = []
    if haskey(project, "authors")
        for author in project["authors"]
            if haskey(author, "orcid")
                person = get_person_from_orcid(author["orcid"])
                push!(authors, person["name_with_email"])
                if haskey(author, "name")
                    person["full_name"] == author["name"] ||
                        @warn "Name mismatch between ORCID and Project.toml: " *
                              "$(person["full_name"]) ≠ $(author["name"])"
                end
                if haskey(author, "email") && haskey(person, "email")
                    person["email"] == author["email"] ||
                        @warn "Email mismatch between ORCID and Project.toml: " *
                              "$(person["email"]) ≠ $(author["email"])"
                end
            elseif haskey(author, "name")
                if haskey(author, "email")
                    push!(authors,
                          author["name"] * " <" * author["email"] * ">")
                else
                    push!(authors, author["name"])
                end
            else
                @warn "Missing name and ORCID in authors block"
            end
        end
    else
        authors = project["author"]
        for author in authors
            name = replace(author, r" *<.*> *$" => "")
            if '<' ∈ author
                email = replace(author, r"^.*<([^>]+)>.*$" => s"\1")
                push!(author_data, Dict("name" => name, "email" => email))
            else
                push!(author_data, Dict("name" => name))
            end
        end
        project["authors"] = author_data
    end

    replace_authors = true
    if !isempty(authors)
        proj_authors = project["author"]
        for author in authors
            if author ∉ proj_authors
                if lowercase(author) ∈ lowercase.(proj_authors)
                    @info "Changing case of $author in Project.toml"
                elseif any(occursin.(author, proj_authors))
                    @info "Completing $author in Project.toml"
                else
                    @error "Author $author not in $proj_authors"
                    replace_authors = false
                end
            end
        end

        if length(authors) < length(proj_authors)
            @error "Author mismatch between $authors and $proj_authors"
            replace_authors = false
        end
    end

    if replace_authors
        project["author"] = authors
    end

    haslicense = false
    license = nothing

    if haskey(project, "license")
        proj_license = project["license"]["SPDX"]
        cm_license = "https://spdx.org/licenses/" * proj_license
        if haskey(codemeta, "license")
            if codemeta["license"] == cm_license
                haslicense = true
                license = proj_license
            else
                @error "License mismatch between Project.toml and codemeta.json: " *
                       "$(codemeta["license"]) ≠ $cm_license"
            end
        else
            codemeta["license"] = cm_license
            haslicense = true
            license = proj_license
        end
    else
        if haskey(codemeta, "license")
            project["license"] = Dict("SPDX" => replace(codemeta["license"],
                                                        "https://spdx.org/licenses/" => ""))
            haslicense = true
            license = project["license"]["SPDX"]
        else
            @warn "No license metadata"
        end
    end

    open_license = nothing
    if haslicense
        url = "https://spdx.org/licenses/$license.json"
        headers = ["Accept" => "application/json"]
        response = HTTP.get(url, headers)

        if response.status == 200
            just_names = replace.(project["author"], r" *<[^>]+> *" => "")
            name_list = join(just_names, ", ", " and ")
            json = JSON.parse(String(response.body))
            open_license = json["isOsiApproved"]
            text = json["licenseText"]
            content = replace(text,
                              r"<year>"i => years,
                              r"<owners?>"i => name_list,
                              r"<copyright holders?>"i => name_list,
                              r"<Owner Organization Name>"i => name_list,
                              r"<Asset Owner>"i => name_list,
                              r"<HOLDERS?>"i => name_list,
                              r"<name of author>"i => name_list,
                              r"<author's name or designee>"i => name_list)
            open("LICENSE", "w") do file
                return write(file, content)
            end
            rm("LICENSE.md", force = true)
        end
    end

    open("Project.toml", "w") do io
        return TOML.print(io, project)
    end
    project = read_project()
    open("Project.toml", "w") do io
        return TOML.print(io, project)
    end

    cm_authors = get(codemeta, "author", OrderedDict{String, Any}[])
    proj_authors = project["authors"]
    cm_from_proj = OrderedDict{String, Any}[]
    for author in proj_authors
        if haskey(author, "orcid")
            person = get_person_from_orcid(author["orcid"])
            dict = OrderedDict{String, Any}("type" => "Person")
            if haskey(person, "givenName")
                dict["givenName"] = person["givenName"]
            end
            if haskey(person, "familyName")
                dict["familyName"] = person["familyName"]
            end
            if haskey(person, "email")
                dict["email"] = person["email"]
            end
            if haskey(person, "orcid")
                dict["id"] = person["pid"]
            end
            if haskey(author, "affiliation")
                dict["affiliation"] = OrderedDict{String, String}[]
                for org in author["affiliation"]
                    d = OrderedDict("type" => "Organization")
                    if haskey(org, "ror")
                        ror = org["ror"]
                        info = get_organisation_from_ror(ror)
                        d["name"] = info["name"]
                        d["identifier"] = info["pid"]
                    else
                        d["name"] = org["name"]
                    end
                    push!(dict["affiliation"], d)
                end
            end
            push!(cm_from_proj, dict)
        end
    end

    if isempty(cm_authors)
        @info "Filling codemeta authors from Project.toml"
        codemeta["author"] = cm_from_proj
    else
        if length(cm_authors) == length(cm_from_proj)
            ids = [dict["id"] for dict in cm_from_proj]
            for dict in cm_authors
                if get(dict, "id", nothing) ∉ ids
                    @error "$dict not found in Project.toml"
                end
            end
        else
            @error "Mismatch between Project.toml and codemeta.json authors: " *
                   "$(cm_authors) ≠ $(proj_authors)"
        end
    end

    if haskey(codemeta, "continuousIntegration")
        if haskey(codemeta, "codemeta:contIntegration")
            if codemeta["continuousIntegration"] ≠
               codemeta["codemeta:contIntegration"]["id"]
                @error "Mismatch between continuousIntegration and codemeta:contIntegration, " *
                       "$(codemeta["continuousIntegration"]) ≠ $(codemeta["codemeta:contIntegration"]["id"])"
            end
        else
            codemeta["codemeta:contIntegration"] = Dict("id" => codemeta["continuousIntegration"])
        end
    else
        if haskey(codemeta, "codemeta:contIntegration")
            codemeta["continuousIntegration"] = codemeta["codemeta:contIntegration"]["id"]
        elseif isfile(".github/workflows/testing.yaml")
            @info "Using .github/workflows/testing.yaml for CI"
            codemeta["continuousIntegration"] = urls[repo_index] *
                                                "/actions/workflows/testing.yaml"
            codemeta["codemeta:contIntegration"] = Dict("id" => codemeta["continuousIntegration"])
        elseif isfile(".github/workflows/CI.yaml")
            @info "Using .github/workflows/CI.yaml for CI"
            codemeta["continuousIntegration"] = urls[repo_index] *
                                                "/actions/workflows/CI.yaml"
            codemeta["codemeta:contIntegration"] = Dict("id" => codemeta["continuousIntegration"])
        elseif isdir(".github/workflows")
            @warn "CI not found in codemeta.json, but .github/workflows exists"
        end
    end

    keywords = sort(get(codemeta, "keywords", ["ecology", "EcoJulia", "julia"]))
    codemeta["keywords"] = keywords

    open("codemeta.json", "w") do io
        return JSON.print(io, codemeta, 4)
    end

    crosswalk = OrderedDict{String, Any}()
    crosswalk["title"] = codemeta["name"]
    if haskey(codemeta, "description")
        crosswalk["description"] = codemeta["description"]
    end
    crosswalk["upload_type"] = "software"
    crosswalk["creators"] = []
    for author in codemeta["author"]
        dict = OrderedDict{String, String}()
        dict["name"] = "$(author["familyName"]), $(author["givenName"])"
        if haskey(author, "id")
            dict["orcid"] = replace(author["id"], "https://orcid.org/" => "")
        end
        if haskey(author, "affiliation")
            affiliations = author["affiliation"]
            affiliation = affiliations isa Vector ? first(affiliations) :
                          affiliations
            if haskey(affiliation, "name")
                dict["affiliation"] = affiliation["name"]
            end
            if haskey(affiliation, "identifier")
                dict["ror"] = replace(affiliation["identifier"],
                                      "https://ror.org/" => "")
            end
        end
        push!(crosswalk["creators"], dict)
    end
    if !isnothing(open_license)
        crosswalk["access_right"] = open_license ? "open" : "closed"
    end
    crosswalk["license"] = project["license"]["SPDX"]
    dict = OrderedDict{String, String}()
    dict["scheme"] = "url"
    dict["identifier"] = codemeta["codeRepository"]
    dict["relation"] = "isOriginalFormOf"
    crosswalk["related_identifiers"] = [dict]
    crosswalk["keywords"] = codemeta["keywords"]

    open(".zenodo.json", "w") do io
        return JSON.print(io, crosswalk, 4)
    end

    # Recursively walk through the directory
    for (root, _, files) in walkdir(".")
        if !startswith(root, "./.git")
            for file in files
                if endswith(file, ".jl")
                    jl_file = joinpath(root, file)
                    data = readlines(jl_file)
                    if startswith(data[1], "# SPDX-License-Identifier:")
                        data[1] = "# SPDX-License-Identifier: $(project["license"]["SPDX"])"
                    else
                        pushfirst!(data, "")
                        pushfirst!(data,
                                   "# SPDX-License-Identifier: $(project["license"]["SPDX"])")
                    end
                    open(jl_file, "w") do io
                        return println.(io, data)
                    end
                end
            end
        end
    end
    return nothing
end

end
