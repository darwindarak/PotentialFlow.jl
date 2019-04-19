using NBInclude

@testset "Jupyter notebooks" begin
    notebooks = ["Pitching Plate - K07",
                 "Point Source Demo",
                 "Translating Plate - 20°",
                 "Translating Plate - 60°",
                 "Vortex Sheet Roll-up",
                 "Moving bodies",
                 "Airfoil in uniform flow"]

    notebook_dir = joinpath(splitdir(@__FILE__())[1], "../binder/notebooks")

    @testset "$notebook" for notebook in notebooks
        @nbinclude(joinpath(notebook_dir, "$notebook.ipynb"); softscope=true)
        @test true
    end
end
