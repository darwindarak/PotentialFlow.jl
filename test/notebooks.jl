using NBInclude

@testset "Jupyter notebooks" begin

    notebook_dir = splitdir(@__FILE__())[1]

    for (root, dirs, files) in walkdir("$notebook_dir/../examples")
        if !contains(root, ".ipynb_checkpoints")
            notebooks = filter(f -> splitext(f)[2] == ".ipynb", files)
            @testset "$notebook" for notebook in notebooks
                nbinclude(joinpath(root, notebook))
                @test true
            end
        end

    end
end
