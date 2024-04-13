"""
    save_epi(file, case_name)

Save the EpiSystem and other inputs generated in files in "test/canonical/" .
These variables are needed to call "simulate!" and are saved as JLSO files

# Arguments
- `file`: file path of a file in "test/canonical/".
- `epi_path`: path to save EpiSystem and other inputs for the simulate run
"""
function save_epi(file, epi_path)
    println("save episystem for case: ", epi_path)
    # prepare saving arguments that will be used in canonical tests through `simulate!`
    #   or `simulate_record!`. These variables are passed into files being included later
    # NOTE: declare `global` variables as `include` will always evaluate files in the global
    #   scope
    global do_save = true
    global save_path = epi_path
    # run file to save epi
    return include(file)
end
