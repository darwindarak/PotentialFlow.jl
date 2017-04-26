using NBInclude

@testset "Jupyter notebooks" begin

    notebook_dir = splitdir(@__FILE__())[1]

    for (root, dirs, files) in walkdir("$notebook_dir/../examples")
        if !contains(root, ".ipynb_checkpoints")
            for file in files
                if splitext(file)[2] == ".ipynb"
                    nbinclude(joinpath(root, file))
                    @test true
                end
            end
        end
    end

end
