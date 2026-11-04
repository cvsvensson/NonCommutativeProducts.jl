# Ranked summary of an AirspeedVelocity run in results/, for the benchmark PR comment:
#     julia benchmark/ci_report.jl <base-rev> <head-rev>
using AirspeedVelocity: load_results

# On recent PRs that leave the code unchanged, no time moved beyond 1.06× and allocated bytes didn't move.
# COARSE keys take tens of ns, where a median has moved from 20 ns to 30 ns.
const TIME_RATIO, COARSE_RATIO, MEMORY_RATIO = 1.15, 1.6, 1.02
const COARSE = ["abstract_eltype/first_segment_range", "abstract_eltype/split_by_group_single", "vector/inner_zero_sum"]
time_ratio(name) = name in COARSE ? COARSE_RATIO : TIME_RATIO

format_time(t) = t < 1.0e-6 ? "$(round(t * 1.0e9, sigdigits = 3)) ns" : t < 1.0e-3 ? "$(round(t * 1.0e6, sigdigits = 3)) µs" :
    t < 1 ? "$(round(t * 1.0e3, sigdigits = 3)) ms" : "$(round(t, sigdigits = 3)) s"
format_bytes(b) = b < 2^10 ? "$(round(Int, b)) B" : b < 2^20 ? "$(round(b / 2^10, sigdigits = 3)) KiB" :
    b < 2^30 ? "$(round(b / 2^20, sigdigits = 3)) MiB" : "$(round(b / 2^30, sigdigits = 3)) GiB"

ratio_beyond(a, b, threshold) = max(a / b, b / a) > threshold

function rows(base, head, key, scale, threshold)
    names = sort!([name for name in intersect(keys(base), keys(head)) if !isnothing(get(base[name], key, nothing))])
    return map(names) do name
        baseline, target = scale * base[name][key], scale * head[name][key]
        (; name, baseline, target, changed = ratio_beyond(target, baseline, threshold(name)), threshold = threshold(name))
    end
end

# Prints a count of the regressed, improved and unchanged rows, then a Markdown table of the changed
# ones: regressions first, each group ranked by ratio.
function print_changes(io, metric, rows, format)
    ratio(r) = r.target / r.baseline
    changed = filter(r -> r.changed, rows)
    regressions = sort!(filter(r -> ratio(r) > 1, changed); by = ratio, rev = true)
    improvements = sort!(filter(r -> ratio(r) < 1, changed); by = ratio)
    println(io, "**$metric**: $(length(regressions)) regressed, $(length(improvements)) improved, $(length(rows) - length(changed)) unchanged\n")
    isempty(changed) && return
    columns = ["ratio" => r -> "$(round(ratio(r), sigdigits = 3))×", "threshold" => r -> "$(r.threshold)×", "benchmark" => r -> "`$(r.name)`", "baseline" => r -> format(r.baseline), "target" => r -> format(r.target)]
    cells = [[first.(columns)]; [[f(r) for f in last.(columns)] for r in [regressions; improvements]]]
    widths = [maximum(c -> textwidth(c[j]), cells) for j in eachindex(columns)]
    line(c) = println(io, "| ", join(rpad.(c, widths), " | "), " |")
    line(cells[1])
    line(["-"^w for w in widths])
    foreach(line, cells[2:end])
    return println(io)
end

base, head = ARGS
results = load_results("NonCommutativeProducts", [base, head]; input_dir = "results")
println("## Benchmarks\n")
print_changes(stdout, "Time", rows(results[base], results[head], "median", 1.0e-9, time_ratio), format_time)
print_changes(stdout, "Memory", rows(results[base], results[head], "memory", 1, _ -> MEMORY_RATIO), format_bytes)
println(
    "A median time counts as changed beyond $(TIME_RATIO)×, or $(COARSE_RATIO)× for the keys of tens of ns listed in `ci_report.jl`; ",
    "memory beyond $(MEMORY_RATIO)×. Ratios are PR / `$base`.\n"
)
