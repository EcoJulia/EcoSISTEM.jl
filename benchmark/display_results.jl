using PkgBenchmark
using Markdown

function display_result(result::BenchmarkResults)
    md = sprint(export_markdown, result)
    md = replace(md, ":x:" => "❌")
    md = replace(md, ":white_check_mark:" => "✅")
    return display(Markdown.parse(md))
end

function get_result(result_file)
    result = PkgBenchmark.readresults(joinpath(@__DIR__, result_file))
    return display_result(result)
end
