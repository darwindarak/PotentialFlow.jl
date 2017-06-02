var documenterSearchIndex = {"docs": [

{
    "location": "index.html#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "index.html#VortexModel-1",
    "page": "Home",
    "title": "VortexModel",
    "category": "section",
    "text": "a scaffolding for building vortex modelsThe main goal of this library is to remove as much boilerplate code as possible from vortex modeling codes. The core operation in vortex models is simulating the dynamics of various interacting vortex elements. In general, the simulation comes down to computing the velocities of the vortex elements then applying some time-marching scheme to evolve the system forward in time. With this in mind, we want to construct a library that makes iteasy to define new vortex types and behaviors\nstraightforward for users to set up a system the vortex elements\nintuitive to probe the state of any vortex element in the system\neasy to define new time-marching schemes to fit the users needs"
},

{
    "location": "index.html#Installation-1",
    "page": "Home",
    "title": "Installation",
    "category": "section",
    "text": "This package requires Julia 0.6- and above. It is not a registered package, so it should be installed with:julia> Pkg.clone(\"git@github.com:darwindarak/VortexModel.jl.git\")Since it is still under heavy development, you should runjulia> Pkg.test(\"VortexModel\") # might take some timeto make sure things are working as intended andjulia> Pkg.update()to get the most recent version of the library and its dependencies.The plots in this documentation are generated using PyPlot.jl. You might want to install that too to follow the examples in the getting started guide or the Jupyter notebooks."
},

{
    "location": "quickstart.html#",
    "page": "Getting Started Guide",
    "title": "Getting Started Guide",
    "category": "page",
    "text": ""
},

{
    "location": "quickstart.html#getting-started-1",
    "page": "Getting Started Guide",
    "title": "Getting Started",
    "category": "section",
    "text": "This getting started guide will introduce the main components of VortexModel.jl. The code examples here should be directly copy-paste-able into the Julia REPL (even with the julia> prompt and sample results).DocTestSetup = quote\n    srand(1)\nend"
},

{
    "location": "quickstart.html#Creating-Vortex-Elements-1",
    "page": "Getting Started Guide",
    "title": "Creating Vortex Elements",
    "category": "section",
    "text": "We start by importing the library and creating a single point vortex with unit circulation located at (1,1):julia> using VortexModel\n\njulia> p = Vortex.Point( 1.0 + 1.0im, 1.0 )\nPoint Vortex: z = 1.0 + 1.0im, Γ = 1.0By convention, the arguments for vortex type constructors are position(s) and circulation(s), followed by any type specific parameters. For example, a vortex blob at the same location as p with a blob radius of 0.1 is created withjulia> Vortex.Blob(1.0 + 1.0im, 1.0, 0.1)\nVortex Blob: z = 1.0 + 1.0im, Γ = 1.0, δ = 0.1We can use Julia's vectorized dot syntax to construct whole arrays of vortex elements. For example, here we create five point vortices and five vortex blobs:julia> N = 5;\n\njulia> zs = Complex.(randn(N), randn(N));\n\njulia> points = Vortex.Point.(zs + 1.5, rand(N))\n5-element Array{VortexModel.Vortex.Points.Point,1}:\n Point Vortex: z = 1.797 + 0.311im, Γ = 0.425\n Point Vortex: z = 1.882 + 2.295im, Γ = 0.773\n Point Vortex: z = 0.902 - 2.267im, Γ = 0.281\n Point Vortex: z = 1.49 + 0.53im, Γ = 0.209\n Point Vortex: z = 0.661 + 0.431im, Γ = 0.251\n\njulia> blobs = Vortex.Blob.(zs - 1.5, rand(N), 0.1)\n5-element Array{VortexModel.Vortex.Blobs.Blob,1}:\n Vortex Blob: z = -1.203 + 0.311im, Γ = 0.02, δ = 0.1\n Vortex Blob: z = -1.118 + 2.295im, Γ = 0.288, δ = 0.1\n Vortex Blob: z = -2.098 - 2.267im, Γ = 0.86, δ = 0.1\n Vortex Blob: z = -1.51 + 0.53im, Γ = 0.077, δ = 0.1\n Vortex Blob: z = -2.339 + 0.431im, Γ = 0.64, δ = 0.1We can mix different vortex types together by grouping them in tuples. For example, a collection of vortex elements consisting of the point vortices and vortex blobs created earlier can be grouped together with:julia> sys = (points, blobs);note: Note\nThe Unicode characters used in the examples can be entered in the Julia REPL (and most text editors with the appropriate plugins) via tab completion..  For example:Γ: \\Gamma<TAB>\nΔ: \\Delta<TAB>\nẋ: x\\dot<TAB>\n🌀: \\:cyclone:<TAB>We can access properties of any vortex element by directly accessing its fields, for example:julia> p.Γ\n1.0However, it is better practice to use accessor methods, such as:julia> Vortex.circulation(p)\n1.0since not all vortex element types store their circulation in a Γ field but all types are required to implement a Vortex.circulation method (also see Vortex.impulse and Vortex.position). These accessor methods, combined with the dot syntax, also make it easier to work with properties of arrays and tuples of vortex elements.julia> Vortex.circulation(points)\n1.939982714228534\n\njulia> Vortex.circulation(blobs)\n1.8849356499471654\n\njulia> Vortex.circulation(sys)\n3.8249183641756996\n\njulia> Vortex.circulation.(blobs)\n5-element Array{Float64,1}:\n 0.0203749\n 0.287702\n 0.859512\n 0.0769509\n 0.640396\n\njulia> Vortex.position.(blobs)\n5-element Array{Complex{Float64},1}:\n -1.20271+0.311111im\n  -1.1176+2.29509im\n -2.09763-2.26709im\n -1.51045+0.529966im\n -2.33903+0.431422im"
},

{
    "location": "quickstart.html#Computing-Vortex-Velocities-1",
    "page": "Getting Started Guide",
    "title": "Computing Vortex Velocities",
    "category": "section",
    "text": "Now that we can create vortex elements, we want to add in some dynamics. The key functions for this are the induce_velocity and induce_velocity! pair and self_induce_velocity!.induce_velocity(target, source) computes the complex velocity that a vortex element(s) source induces on a target. The target can bea complex position\njulia> induce_velocity(0.0 + 0.0im , points)\n0.05610938572529216 - 0.1319030126670981im\n\njulia> induce_velocity(0.0 + 0.0im , sys)\n0.05066830110387291 - 0.04224547600656549im\na vortex element\njulia> induce_velocity(p, sys)\n-0.095439940976663 - 0.024542142467999073im\nan array/tuple of vortex elements\njulia> induce_velocity(points, blobs)\n5-element Array{Complex{Float64},1}:\n -0.00789749+0.0645051im\n  -0.0278927+0.0538741im\n   0.0271037+0.0706032im\n  -0.0111193+0.0675933im\n  -0.0117893+0.078857im\n\njulia> induce_velocity(blobs, sys)\n5-element Array{Complex{Float64},1}:\n  0.0126862+0.0352193im\n  -0.111207-0.0472771im\n  0.0873796-0.0535197im\n -0.0375196+0.031068im\n -0.0279267-0.103821imThe in-place version, induce_velocity!(velocities, targets, source), computes the velocity and writes the results into a pre-allocated data structure. For example:julia> vel_points = zeros(Complex128, length(points))\n5-element Array{Complex{Float64},1}:\n 0.0+0.0im\n 0.0+0.0im\n 0.0+0.0im\n 0.0+0.0im\n 0.0+0.0im\n\njulia> induce_velocity!(vel_points, points, blobs);\n\njulia> vel_points\n5-element Array{Complex{Float64},1}:\n -0.00789749+0.0645051im\n  -0.0278927+0.0538741im\n   0.0271037+0.0706032im\n  -0.0111193+0.0675933im\n  -0.0117893+0.078857imTo make it easier to allocate velocities for more complex collections of vortex elements, the library provides the allocate_velocity function:julia> vels = allocate_velocity(sys);\n\njulia> typeof(vels)\nTuple{Array{Complex{Float64},1},Array{Complex{Float64},1}}The code above created a tuple containing two arrays of velocities, corresponding to the structure of sys. Similarly, there is also the reset_velocity!(velocities, sources) function, which resizes the entries in velocities to match the structure of sources if necessary, then sets all velocities to zero. We can compute the velocity that a source induces on the entire points/blobs system with:julia> src = Vortex.Point(1.0, 1.0);\n\njulia> induce_velocity!(vels, sys, src)\n(Complex{Float64}[-0.067601+0.173242im, -0.0604154+0.023228im, 0.0700725-0.00301774im, -0.162041+0.149685im, -0.228068-0.179224im], Complex{Float64}[-0.0100056-0.0708409im, -0.0374576-0.0345609im, 0.0244871-0.033458im, -0.0128124-0.0606923im, -0.00605748-0.0468824im])If we want the velocity that the points/blobs system induces on itself, we can callreset_velocity!(vels, sys)\ninduce_velocity!(vels[1], points, points)\ninduce_velocity!(vels[1], points, src)\ninduce_velocity!(vels[2], blobs, src)\ninduce_velocity!(vels[2], blobs, blobs)This becomes difficult to keep track of when sys gets larger or more complicated (e.g. nested collection of elements). Instead, we can use the self_induce_velocity! function, which takes care of applying all the pairwise interactions (recursively if need be):julia> reset_velocity!(vels, sys);\n\njulia> self_induce_velocity!(vels, sys);"
},

{
    "location": "quickstart.html#Time-Marching-1",
    "page": "Getting Started Guide",
    "title": "Time Marching",
    "category": "section",
    "text": "using VortexModel\nusing Gadfly\nimport Colors: colormap, alphacolor\nsrand(1)\n\nfunction plot_system(sys)\n    plot(x = real.(Vortex.position.(vcat(sys...))),\n         y = imag.(Vortex.position.(vcat(sys...))),\n         color = Vortex.circulation.(vcat(sys...)),\n         Coord.cartesian(fixed=true),\n         Guide.colorkey(\"Γ\"),\n         Scale.color_continuous(colormap=Scale.lab_gradient(colormap(\"reds\")...)),\n         style(grid_line_width=0mm, highlight_width=0mm))\nendNow that we compute the velocities of a system of vortex elements, we can march the system forward in time to simulate its behavior. As an example, we will simulate of two clusters of vortex blobs merging.N = 200\nzs = Complex.(0.5randn(N), 0.5randn(N))\nΓs  = @. exp(-4abs2(zs))\ncluster₁ = Vortex.Blob.(zs + 1, Γs, 0.01)\ncluster₂ = Vortex.Blob.(zs - 1, Γs, 0.01)\n\nsys = (cluster₁, cluster₂)\nvels = allocate_velocity(sys)\nplot_system(sys)\ndraw(SVGJS(\"initial_clusters.svg\", 6inch, 4inch), ans); nothing # hidewarning: Warning\nFunctions for plotting vortex elements are still waiting for a couple more issues to be fixed on Plots.jl.  For now, we can use Gadfly.jl directly as follows:using Gadfly\nplot(x = real.(Vortex.position.(vcat(sys...))),\n     y = imag.(Vortex.position.(vcat(sys...))),\n     color = Vortex.circulation.(vcat(sys...)),\n     Coord.cartesian(fixed=true),\n     Guide.colorkey(\"Γ\"),\n     Scale.color_continuous(colormap=Scale.lab_gradient(colormap(\"reds\")...)),\n     style(grid_line_width=0mm, highlight_width=0mm))Alternatively, we can use PyPlot.jl with something like:using PyPlot\nfor cluster in sys\n    scatter(real.(Vortex.position.(cluster)),\n            imag.(Vortex.position.(cluster)),\n            c = Vortex.circulation.(cluster),\n            vmin = 0, vmax = 1, alpha = 0.7,\n            cmap = PyPlot.get_cmap(\"Reds\"))\nend\ncolorbar()\naxis(:scaled)\naxis([-3,3,-3,3])<object data=\"initial_clusters.svg\" type=\"image/svg+xml\"></object>Given an array or tuple of vortex elements and their velocities, we can compute their positions after some time interval with the advect!(x₊, x, ẋ, Δt) function, wherex₊ is where the new states are stored\nx is the current state\nΔt is the time interval\nẋ is the velocity.In our case, we will let x₊ and x both be set to sys:Δt = 0.01\nfor t in 0:Δt:1.0\n    reset_velocity!(vels, sys)\n    self_induce_velocity!(vels, sys)\n    advect!(sys, sys, vels, Δt)\nend\nplot_system(sys)\ndraw(SVGJS(\"final_clusters.svg\", 6inch, 4inch), ans); nothing # hide<object data=\"final_clusters.svg\" type=\"image/svg+xml\"></object>"
},

{
    "location": "elements.html#",
    "page": "Vortex Elements",
    "title": "Vortex Elements",
    "category": "page",
    "text": ""
},

{
    "location": "elements.html#Vortex-Elements-1",
    "page": "Vortex Elements",
    "title": "Vortex Elements",
    "category": "section",
    "text": "DocTestSetup = quote\nusing VortexModel\nsrand(1)\nendThe library currently has four built-in vortex types:Vortex.Point\nVortex.Blob\nVortex.Sheet\nVortex.Plate (at the moment, there can only be one plate in the fluid at at time)Most functions in the library that act on vortex elements can take either a single vortex element, or a collection of elements. These collections can be represented as an array or a tuple. Arrays should be used when the elements are the same type, for example:julia> points = Vortex.Point.(rand(Complex128, 5), rand(5))\n5-element Array{VortexModel.Vortex.Points.Point,1}:\n Point Vortex: z = 0.236 + 0.347im, Γ = 0.556\n Point Vortex: z = 0.313 + 0.008im, Γ = 0.437\n Point Vortex: z = 0.489 + 0.211im, Γ = 0.425\n Point Vortex: z = 0.952 + 1.0im, Γ = 0.773\n Point Vortex: z = 0.252 + 0.987im, Γ = 0.281\n\njulia> Vortex.impulse(points)\n1.3362266530178137 - 1.2821936908564113im\n\njulia> blobs = [Vortex.Blob(rand(Complex128), rand(), 0.1) for i in 1:5]\n5-element Array{VortexModel.Vortex.Blobs.Blob,1}:\n Vortex Blob: z = 0.209 + 0.251im, Γ = 0.02, δ = 0.1\n Vortex Blob: z = 0.288 + 0.86im, Γ = 0.077, δ = 0.1\n Vortex Blob: z = 0.64 + 0.874im, Γ = 0.279, δ = 0.1\n Vortex Blob: z = 0.751 + 0.645im, Γ = 0.078, δ = 0.1\n Vortex Blob: z = 0.848 + 0.086im, Γ = 0.553, δ = 0.1\n\njulia> Vortex.impulse(blobs)\n0.41217890550975256 - 0.7325028967929701imKnowing that every element has the same type allows the compiler to perform more aggressive optimizations. Tuples are used when we want to mix and match different vortex types. For example:julia> sys = (points, blobs);\n\njulia> Vortex.impulse(sys)\n1.7484055585275664 - 2.0146965876493814imThis rest of this page documents the data types that represent these elements and some key functions that act on them. For more detailed examples, please refer to the Jupyter notebooks."
},

{
    "location": "elements.html#VortexModel.Vortex.Points.Point",
    "page": "Vortex Elements",
    "title": "VortexModel.Vortex.Points.Point",
    "category": "Type",
    "text": "Vortex.Point <: Vortex.PointSource\n\nAn immutable structure representing a point vortex\n\nFields\n\nz: position\nΓ: circulation\n\n\n\n"
},

{
    "location": "elements.html#VortexModel.Vortex.Blobs.Blob",
    "page": "Vortex Elements",
    "title": "VortexModel.Vortex.Blobs.Blob",
    "category": "Type",
    "text": "Vortex.Blob <: Vortex.PointSource\n\nAn immutable structure representing a vortex blob\n\nFields\n\nz: position\nΓ: circulation\nδ: blob radius\n\n\n\n"
},

{
    "location": "elements.html#VortexModel.Vortex.Sheets.Sheet",
    "page": "Vortex Elements",
    "title": "VortexModel.Vortex.Sheets.Sheet",
    "category": "Type",
    "text": "Vortex.Sheet <: Vortex.CompositeSource\n\nA vortex sheet represented by vortex blob control points\n\nFields\n\nblobs: the underlying array of vortex blobs\nΓs: the cumulated sum of circulation starting from the first control point\nδ: the blob radius of all the vortex blobs\nzs: a mapped array that accesses the position of each control point\n\nConstructors:\n\nSheet(blobs, Γs, δ)\nSheet(zs, Γs, δ) where zs is an array of positions for the control points\n\n\n\n"
},

{
    "location": "elements.html#VortexModel.Vortex.Plates.Plate",
    "page": "Vortex Elements",
    "title": "VortexModel.Vortex.Plates.Plate",
    "category": "Type",
    "text": "Vortex.Plate <: VortexCompositeSource\n\nAn infinitely thin, flat plate, represented as a bound vortex sheet\n\nFields\n\nL\nchord length\nc\ncentroid\nα\ncentroid velocity\nΓ\ntotal circulation\nN\nnumber of control points\nss\nnormalized positions (within [-1, 1]) of the control points\nzs\ncontrol point coordinates\nA\nChebyshev coefficients of the normal component of velocity induced along the plate by ambient vorticity\nC\nChebyshev coefficients of the velocity induced along the plate by ambient vorticity\nB₀\nzeroth Chebyshev coefficient associated with body motion\nB₁\nfirst Chebyshev coefficient associated with body motion\ndchebt!\nPreplanned discrete Chebyshev transform\n\nConstructors\n\nPlate(N, L, c, α)\n\n\n\n"
},

{
    "location": "elements.html#Built-in-Vortex-Types-1",
    "page": "Vortex Elements",
    "title": "Built-in Vortex Types",
    "category": "section",
    "text": "Vortex.Point\nVortex.Blob\nVortex.Sheet\nVortex.Plate"
},

{
    "location": "elements.html#VortexModel.Vortex.position",
    "page": "Vortex Elements",
    "title": "VortexModel.Vortex.position",
    "category": "Function",
    "text": "Vortex.position(src::PointSource)\n\nReturns the complex position of a PointSource type vortex element This is a required method for all subtypes of PointSource.\n\nExample\n\njulia> point = Vortex.Point(1.0 + 0.0im, 1.0);\n\njulia> Vortex.position(point)\n1.0 + 0.0im\n\njulia> points = Vortex.Point.([1.0im, 2.0im], 1.0);\n\njulia> Vortex.position.(points)\n2-element Array{Complex{Float64},1}:\n 0.0+1.0im\n 0.0+2.0im\n\n\n\n"
},

{
    "location": "elements.html#VortexModel.Vortex.circulation",
    "page": "Vortex Elements",
    "title": "VortexModel.Vortex.circulation",
    "category": "Function",
    "text": "Vortex.circulation(src)\n\nReturns the total circulation contained in src This is a required method for all vortex types.\n\nExample\n\njulia> points = Vortex.Point.([1.0im, 2.0im], [1.0, 2.0]);\n\njulia> Vortex.circulation(points[1])\n1.0\n\njulia> Vortex.circulation(points)\n3.0\n\njulia> Vortex.circulation.(points)\n2-element Array{Float64,1}:\n 1.0\n 2.0\n\n\n\n"
},

{
    "location": "elements.html#VortexModel.Vortex.impulse",
    "page": "Vortex Elements",
    "title": "VortexModel.Vortex.impulse",
    "category": "Function",
    "text": "Vortex.impulse(src)\n\nReturn the aerodynamic impulse of src about (0,0):\n\nP = int boldsymbolx times boldsymbolomegamathrmdA\n\nThis is a required method for all vortex types.\n\nExample\n\njulia> sys = (Vortex.Point(1.0im, π), Vortex.Blob(1.0im, -π, 0.1));\n\njulia> Vortex.impulse(sys[1])\n3.141592653589793 + 0.0im\n\njulia> Vortex.impulse(sys)\n0.0 + 0.0im\n\n\n\n"
},

{
    "location": "elements.html#VortexModel.Vortex.advect",
    "page": "Vortex Elements",
    "title": "VortexModel.Vortex.advect",
    "category": "Function",
    "text": "advect(src::PointSource, velocity::Complex128, Δt)\n\nReturn a new vortex element that represents src advected by velocity over Δt If this method is implemented by any type T <: PointSource, then an array of type AbstractArray{T} can be passed in the first two arguments of advect!.\n\nExample\n\njulia> point = Vortex.Point(1.0 + 0.0, 1.0);\n\njulia> advect(point, 1.0im, 1e-2)\nPoint Vortex: z = 1.0 + 0.01im, Γ = 1.0\n\n\n\n"
},

{
    "location": "elements.html#VortexModel.Vortex.advect!",
    "page": "Vortex Elements",
    "title": "VortexModel.Vortex.advect!",
    "category": "Function",
    "text": "advect!(srcs₊, srcs₋, vels, Δt)\n\nMoves the elements in srcs₋ by their corresponding velocity in vels over the interval Δt and store the results in src₊ srcs₋ and srcs₊ can be either a array of vortex elements or a tuple.\n\nExample\n\njulia> points₋ = [Vortex.Point(x + 0im, 1.0) for x in 1:5];\n\njulia> points₊ = Vector{Vortex.Point}(5);\n\njulia> vels = [ y*im for y in 1.0:5 ];\n\njulia> advect!(points₊, points₋, vels, 1e-2)\n\njulia> points₊\n5-element Array{VortexModel.Vortex.Points.Point,1}:\n Point Vortex: z = 1.0 + 0.01im, Γ = 1.0\n Point Vortex: z = 2.0 + 0.02im, Γ = 1.0\n Point Vortex: z = 3.0 + 0.03im, Γ = 1.0\n Point Vortex: z = 4.0 + 0.04im, Γ = 1.0\n Point Vortex: z = 5.0 + 0.05im, Γ = 1.0\n\n\n\n"
},

{
    "location": "elements.html#Vortex-Properties-1",
    "page": "Vortex Elements",
    "title": "Vortex Properties",
    "category": "section",
    "text": "Vortex.position\nVortex.circulation\nVortex.impulse\nadvect\nadvect!"
},

{
    "location": "elements.html#VortexModel.Vortex.Sheets.append_segment!",
    "page": "Vortex Elements",
    "title": "VortexModel.Vortex.Sheets.append_segment!",
    "category": "Function",
    "text": "Vortex.Sheets.append_segment!(sheet::Sheet, z, Γ)\n\nAppend a new segment with circulation Γ extending from the end of the sheet to z.\n\nExample\n\njulia> sheet = Vortex.Sheet(0:0.1:1, 0.0:10, 0.2)\nVortex Sheet: L ≈ 1.0, Γ = 10.0, δ = 0.2\n\njulia> sheet.blobs[end]\nVortex Blob: z = 1.0 + 0.0im, Γ = 0.5, δ = 0.2\n\njulia> Vortex.Sheets.append_segment!(sheet, 1.1, 2.0)\n\njulia> sheet\nVortex Sheet: L ≈ 1.1, Γ = 12.0, δ = 0.2\n\njulia> sheet.blobs[end]\nVortex Blob: z = 1.1 + 0.0im, Γ = 1.0, δ = 0.2\n\n\n\n"
},

{
    "location": "elements.html#VortexModel.Vortex.Sheets.truncate!",
    "page": "Vortex Elements",
    "title": "VortexModel.Vortex.Sheets.truncate!",
    "category": "Function",
    "text": "Vortex.Sheets.truncate!(sheet, n::Int)\n\nRemove segments 0:n from sheet, and return the circulation in those segments.\n\nExample\n\njulia> sheet = Vortex.Sheet(0:0.1:1, 0.0:10, 0.2)\nVortex Sheet: L ≈ 1.0, Γ = 10.0, δ = 0.2\n\njulia> Vortex.Sheets.truncate!(sheet, 5)\n4.0\n\n\n\n"
},

{
    "location": "elements.html#VortexModel.Vortex.Sheets.redistribute_points!",
    "page": "Vortex Elements",
    "title": "VortexModel.Vortex.Sheets.redistribute_points!",
    "category": "Function",
    "text": "Vortex.Sheets.redistribute_points!(sheet, zs, Γs)\n\nReturns the modified sheet with replacement control points at positions zs and strength Γs.\n\njulia> sheet = Vortex.Sheet(0:0.1:1, 0.0:10, 0.2)\nVortex Sheet: L ≈ 1.0, Γ = 10.0, δ = 0.2\n\njulia> sys = (sheet,)\n(Vortex Sheet: L ≈ 1.0, Γ = 10.0, δ = 0.2,)\n\njulia> Vortex.Sheets.redistribute_points!(sheet, 0:0.2:2, 0.0:0.5:5)\nVortex Sheet: L ≈ 2.0, Γ = 5.0, δ = 0.2\n\njulia> sys[1]\nVortex Sheet: L ≈ 2.0, Γ = 5.0, δ = 0.2\n\n\n\n"
},

{
    "location": "elements.html#VortexModel.Vortex.Sheets.remesh",
    "page": "Vortex Elements",
    "title": "VortexModel.Vortex.Sheets.remesh",
    "category": "Function",
    "text": "Vortex.Sheets.remesh(sheet, Δs::Float64 , params::Tuple = ())\n\nUniformly redistribute the control points of the sheet to have a nominal spacing of Δs. Material quantities that should be redistributed along with the control points can be passed in as elements of params.\n\nReturns the tuple (z₌, Γ₌, L [, p₌]) where\n\nz₌ is an array with the positions of the uniformly distributed points\nΓ₌ is circulation interpolated onto z₌\nL is total length of the sheet\np₌ is a tuple containing the material quantities from params interpolated onto z₌\n\nExample\n\njulia> sheet = Vortex.Sheet(0:0.1:1, 0.0:10, 0.2)\nVortex Sheet: L ≈ 1.0, Γ = 10.0, δ = 0.2\n\njulia> age = collect(10.0:-1:0);\n\njulia> Vortex.Sheets.remesh(sheet, 0.2, (age, ))\n(Complex{Float64}[0.0+0.0im, 0.25+0.0im, 0.5+0.0im, 0.75+0.0im, 1.0+0.0im], [0.0, 2.5, 5.0, 7.5, 10.0], 1.0, ([10.0, 7.5, 5.0, 2.5, 0.0],))\n\n\n\n"
},

{
    "location": "elements.html#VortexModel.Vortex.Sheets.remesh!",
    "page": "Vortex Elements",
    "title": "VortexModel.Vortex.Sheets.remesh!",
    "category": "Function",
    "text": "Vortex.Sheets.remesh!(sheet::Sheet, Δs::Float64, params::Tuple = ())\n\nSame as Vortex.Sheets.remesh, except sheet is replaced internally by a uniformly interpolated control points. Returns the tuple (sheet, L, p₌) where\n\nsheet is the modified sheet\nL is total length of the sheet\np₌ is a tuple containing the material quantities from params interpolated onto the new control points of sheet\n\njulia> sheet = Vortex.Sheet(0:0.1:1, 0.0:10, 0.2)\nVortex Sheet: L ≈ 1.0, Γ = 10.0, δ = 0.2\n\njulia> age = collect(10.0:-1:0);\n\njulia> Vortex.Sheets.remesh!(sheet, 0.2, (age,));\n\njulia> Vortex.position.(sheet.blobs)\n5-element Array{Complex{Float64},1}:\n  0.0+0.0im\n 0.25+0.0im\n  0.5+0.0im\n 0.75+0.0im\n  1.0+0.0im\n\njulia> age\n5-element Array{Float64,1}:\n 10.0\n  7.5\n  5.0\n  2.5\n  0.0\n\n\n\n"
},

{
    "location": "elements.html#VortexModel.Vortex.Sheets.split!",
    "page": "Vortex Elements",
    "title": "VortexModel.Vortex.Sheets.split!",
    "category": "Function",
    "text": "Vortex.Sheets.split!(sheet, n::Int)\n\nRemove segments 0:n from sheet, and return those segments as a new sheet.\n\nExample\n\njulia> sheet = Vortex.Sheet(0:0.1:1, 0.0:10, 0.2)\nVortex Sheet: L ≈ 1.0, Γ = 10.0, δ = 0.2\n\njulia> sheet₋ = Vortex.Sheets.split!(sheet, 5)\nVortex Sheet: L ≈ 0.4, Γ = 4.0, δ = 0.2\n\njulia> sheet\nVortex Sheet: L ≈ 0.6, Γ = 6.0, δ = 0.2\n\n\n\n"
},

{
    "location": "elements.html#VortexModel.Vortex.Sheets.filter!",
    "page": "Vortex Elements",
    "title": "VortexModel.Vortex.Sheets.filter!",
    "category": "Function",
    "text": "Vortex.Sheets.filter!(sheet, Δs, Δf[, params])\n\nRedistribute and filter the control points of a vortex sheet \n\nArguments\n\nsheet: the vortex sheet to be modified\nΔs: the nominal spacing between the uniform points\nΔf: the minimum length scale that the filter should allow to pass through\nparams: an optional tuple of vectors containing material properties\n\nReturns\n\nIf params is passed in, then its vectors will be overwritten by their interpolated values on the new control points, and the function returns the tuple (sheet, params). Otherwise, it returns (sheet, ())\n\n\n\n"
},

{
    "location": "elements.html#VortexModel.Vortex.Sheets.filter_position!",
    "page": "Vortex Elements",
    "title": "VortexModel.Vortex.Sheets.filter_position!",
    "category": "Function",
    "text": "filter_position!(s, Δf, L = arclength(z₌))\n\nFilter out any length scales in s that is smaller than Δf, storing the result back in s. s can be either a vector of complex positions, or a Vortex.Sheet.\n\n\n\n"
},

{
    "location": "elements.html#VortexModel.Vortex.Sheets.arclength",
    "page": "Vortex Elements",
    "title": "VortexModel.Vortex.Sheets.arclength",
    "category": "Function",
    "text": "arclength(s)\n\nCompute the polygonal arc length of s, where s can be either an vector of complex numbers or a Vortex.Sheet.\n\nExample\n\n```jldoctest julia> sheet = Vortex.Sheet(0:0.1:1, 0.0:10, 0.2) Vortex Sheet: L ≈ 1.0, Γ = 10.0, δ = 0.2\n\njulia> Vortex.Sheets.arclength(sheet) 1.0\n\n\n\n"
},

{
    "location": "elements.html#VortexModel.Vortex.Sheets.arclengths",
    "page": "Vortex Elements",
    "title": "VortexModel.Vortex.Sheets.arclengths",
    "category": "Function",
    "text": "arclengths(s)\n\nCumulative sum of the polygonal arc length of s, where s can be either an vector of complex numbers or a Vortex.Sheet.\n\nExample\n\njulia> sheet = Vortex.Sheet(0:0.1:1, 0.0:10, 0.2)\nVortex Sheet: L ≈ 1.0, Γ = 10.0, δ = 0.2\n\njulia> Vortex.Sheets.arclengths(sheet)\n11-element Array{Float64,1}:\n 0.0\n 0.1\n 0.2\n 0.3\n 0.4\n 0.5\n 0.6\n 0.7\n 0.8\n 0.9\n 1.0\n\n\n\n"
},

{
    "location": "elements.html#Methods-on-Vortex-Sheets-1",
    "page": "Vortex Elements",
    "title": "Methods on Vortex Sheets",
    "category": "section",
    "text": "Vortex.Sheets.append_segment!\nVortex.Sheets.truncate!\nVortex.Sheets.redistribute_points!\nVortex.Sheets.remesh\nVortex.Sheets.remesh!\nVortex.Sheets.split!\nVortex.Sheets.filter!\nVortex.Sheets.filter_position!\nVortex.Sheets.arclength\nVortex.Sheets.arclengths"
},

{
    "location": "elements.html#VortexModel.Vortex.Plates.enforce_no_flow_through!",
    "page": "Vortex Elements",
    "title": "VortexModel.Vortex.Plates.enforce_no_flow_through!",
    "category": "Function",
    "text": "enforce_no_flow_through!(p::Plate, motion, elements)\n\nUpdate the plate, p, to enforce the no-flow-through condition given ambient vortex elements, elements, and while moving with kinematics specified by motion.\n\nExample\n\njulia> plate = Vortex.Plate(128, 2.0, 0.0, π/3)\nPlate: N = 128, L = 2.0, c = 0.0 + 0.0im, α = 60.0ᵒ\n       LESP = 0.0, TESP = 0.0\n\njulia> motion = allocate_velocity(plate); motion.ċ = 1.0;\n\njulia> point = Vortex.Point(0.0 + 2im, 1.0);\n\njulia> Vortex.enforce_no_flow_through!(plate, motion, point)\n\njulia> plate\nPlate: N = 128, L = 2.0, c = 0.0 + 0.0im, α = 60.0ᵒ\n       LESP = 1.27, TESP = -1.93\n\n\n\n"
},

{
    "location": "elements.html#VortexModel.Vortex.Plates.vorticity_flux",
    "page": "Vortex Elements",
    "title": "VortexModel.Vortex.Plates.vorticity_flux",
    "category": "Function",
    "text": "vorticity_flux(p::Plate, v₁, v₂,\n               lesp = 0.0, tesp = 0.0,\n               ∂C₁ = Vector{Complex128}(plate.N),\n               ∂C₂ = Vector{Complex128}(plate.N))\n\nReturn strengths of new vortex elements that satisfies edge suction parameters. For a given edge, if the current suction parameter is less than the criticial suction parameter, then no vorticity is released.  If it is higher, however, vorticity will be released so that the suction parameter equals the critical value.\n\nArguments\n\np: the plate\nv₁, v₂: the vortex elements (with unit circulation) that the vorticity flux is going into\nlesp, tesp: the critical leading and trailing edge suction parameters we want to enforce.  By default, both parameters are set to 0.0 to enforce the Kutta condition on both edges.  We can disable vortex shedding from an edge by setting the its critical suction parameter to Inf\n\nReturns\n\nΓ₁, Γ₂: the strengths that the vortex element should have in order to satisfy the edge suction parameters\n∂C₁, ∂C₂: Chebyshev coefficients of the normal velocity induced by the vortex elements Instead of running enforce_bc! with the new vortex elements, we can use this matrix to directly update the Chebyshev coefficients associated with the bound vortex sheet without recomputing all the velocities.\n\nExample\n\nEnforcing the trailing edge Kutta condition with an point vortex at negative infinity:\n\njulia> plate = Vortex.Plate(128, 2.0, 0.0, π/6)\nPlate: N = 128, L = 2.0, c = 0.0 + 0.0im, α = 30.0ᵒ\n       LESP = 0.0, TESP = 0.0\n\njulia> motion = allocate_velocity(plate);\n\njulia> motion.ċ = 1.0;\n\njulia> Vortex.enforce_no_flow_through!(plate, motion, ())\n\njulia> point = Vortex.Point(-Inf, 1.0);\n\njulia> _, Γ, _, _ = Vortex.vorticity_flux(plate, (), point,  Inf);\n\njulia> Γ # should equal -πULsin(α) = -π\n-3.1415926535897927\n\n\n\n"
},

{
    "location": "elements.html#VortexModel.Vortex.Plates.vorticity_flux!",
    "page": "Vortex Elements",
    "title": "VortexModel.Vortex.Plates.vorticity_flux!",
    "category": "Function",
    "text": "vorticity_flux!(p::Plate, v₁, v₂,\n                lesp = 0.0, tesp = 0.0,\n                ∂C₁ = Vector{Complex128}(plate.N),\n                ∂C₂ = Vector{Complex128}(plate.N))\n\nIn-place version of vorticity_flux, except instead of just returning the possible changes in plate Chebyshev coefficients, we modify plate.C with those changes so that no-flow-through is enforced in the presence of v₁ and v₂ with strengths that satisfy the suction parameters.\n\n\n\n"
},

{
    "location": "elements.html#VortexModel.Vortex.Plates.bound_circulation",
    "page": "Vortex Elements",
    "title": "VortexModel.Vortex.Plates.bound_circulation",
    "category": "Function",
    "text": "bound_circulation(plate[, s])\n\nCompute the bound circulation between the trailing edge of the plate to s.\n\ns can be either a single normalized arc length coordinate (between -1 and 1), or a whole array of coordinates.\n\n\n\n"
},

{
    "location": "elements.html#VortexModel.Vortex.Plates.bound_circulation!",
    "page": "Vortex Elements",
    "title": "VortexModel.Vortex.Plates.bound_circulation!",
    "category": "Function",
    "text": "bound_circulation!(Γs, plate[, ss])\n\nCompute the bound circulation between the trailing edge of the plate to ss, then store it in Γs.\n\nIf an array, ss, with normalized arc length coordinates is omitted, then the circulation will be computed at the plate's Chebyshev nodes.\n\n\n\n"
},

{
    "location": "elements.html#VortexModel.Vortex.Plates.surface_pressure",
    "page": "Vortex Elements",
    "title": "VortexModel.Vortex.Plates.surface_pressure",
    "category": "Function",
    "text": "surface_pressure(plate, motion, te_sys, Γs₋, Δt)\n\nCompute the pressure difference across the plate along Chebyshev nodes.\n\nnote: Note\nThe pressure difference across the bound vortex sheet is given by:    p_-^+\n  = -rho left frac12(boldsymbolv^+ + boldsymbolv^-)\n               - boldsymbolv_b\n         right\n         cdot ( boldsymbolgamma cross boldsymbolhatn)\n    +rho fracmathrmdGammamathrmdtwhere rho is the fluid density, boldsymbolv^pm is the velocity on either side of the plate, boldsymbolv_b is the local velocity of the plate, boldsymbolgamma is the bound vortex sheet strength, and Gamma is the integrated circulation. We will compute fracmathrmdGammamathrmdt using finite differences.  So we will need the circulation along the plate from a previous time-step in order to compute the current pressure distribution.  We assume that value of circulation at the trailing edge of the plate is equal the the net circulation of all the vorticity that has been shed from the trailing edge.\n\nArguments\n\nplate: we assume that the Plate structure that is passed in already enforces the no-flow-through condition\nmotion: the motion of the plate used to compute boldsymbolv_b\nte_sys: the system of vortex elements representing the vorticity shed from the trailing edge of the plate\nΓs₋: the circulation along the plate's Chebyshev nodes, this should be equivalent to calling Vortex.circulation(te_sys) .+ Vortex.bound_circulation(plate) from a previous time-step.\nΔt: time-step used to compute ``\\frac{\\mathrm{d}\\Gamma}{\\mathrm{d}t} using finite differences\n\nReturns\n\nΔp: the pressure difference across the plate along Chebyshev nodes\nΓs₊: the circulation along the plate at the current time-step (this value is used in computing the current Δp and can be used as the Γs₋ for computing pressure differences at the next time-step)\n\n\n\n"
},

{
    "location": "elements.html#Methods-on-Plates-1",
    "page": "Vortex Elements",
    "title": "Methods on Plates",
    "category": "section",
    "text": "Vortex.Plates.enforce_no_flow_through!\nVortex.Plates.vorticity_flux\nVortex.Plates.vorticity_flux!\nVortex.Plates.bound_circulation\nVortex.Plates.bound_circulation!\nVortex.Plates.surface_pressure"
},

{
    "location": "elements.html#Index-1",
    "page": "Vortex Elements",
    "title": "Index",
    "category": "section",
    "text": "Pages = [\"elements.md\"]"
},

{
    "location": "velocities.html#",
    "page": "Computing Velocities",
    "title": "Computing Velocities",
    "category": "page",
    "text": ""
},

{
    "location": "velocities.html#Computing-Velocities-1",
    "page": "Computing Velocities",
    "title": "Computing Velocities",
    "category": "section",
    "text": "DocTestSetup = quote\nusing VortexModel\nsrand(1)\nend"
},

{
    "location": "velocities.html#Sources-and-Targets-1",
    "page": "Computing Velocities",
    "title": "Sources and Targets",
    "category": "section",
    "text": "Velocity computations in vortex models essentially boils down to pairwise interactions between sources and targets. We may be interested in how a system of vortex elements induces velocity on at point, at multiple points, on other vortex elements, or on itself.The three key functions for computing velocities areinduce_velocity(target, source)\ninduce_velocity!(velocity, target, source)\nself_induce_velocity!(velocity, source)The ! suffix in the last two function signatures indicate that the velocity argument will be overwritten by the results of the computation.Sources of velocity can be any one of:a single vortex element, e.g.\njulia> src = Vortex.Point(im, 1.0);\n\njulia> induce_velocity(0.0 + 0.0im, src)\n0.15915494309189535 - 0.0im\nan array of homogenous vortex types, e.g.\njulia> srcs = Vortex.Point.([im, 1.0], 1.0);\n\njulia> induce_velocity(0.0 + 0.0im, srcs)\n0.15915494309189535 - 0.15915494309189535im\na tuple of different vortex types, e.g.\njulia> srcs₂ = Vortex.Point.([2im, 2.0], -2.0);\n\njulia> sys = (srcs, srcs₂);\n\njulia> induce_velocity(0.0 + 0.0im, sys)\n0.0 + 0.0imIn the examples above, the target was just complex number 0.0 + 0.0im. However we can also havean array of complex numbers, e.g.\njulia> targets = Complex128.(1:3);\n\njulia> induce_velocity(targets, src)\n3-element Array{Complex{Float64},1}:\n 0.0795775+0.0795775im\n  0.031831+0.063662im\n 0.0159155+0.0477465im\nan array of vortex elements, e.g.\njulia> targets₂ = Vortex.Point.(im*(1.0:3), 1.0);\n\njulia> induce_velocity(targets₂, src)\n3-element Array{Complex{Float64},1}:\n        0.0+0.0im\n  -0.159155+0.0im\n -0.0795775+0.0im\na tuple with any of the above, e.g.\njulia> targets₃ = Vortex.Point.(-3.0:-1, -1.0);\n\njulia> sys = (targets, (targets₂, targets₃));\n\njulia> induce_velocity(sys, src)\n(Complex{Float64}[0.0795775+0.0795775im, 0.031831+0.063662im, 0.0159155+0.0477465im], (Complex{Float64}[0.0+0.0im, -0.159155+0.0im, -0.0795775+0.0im], Complex{Float64}[0.0159155-0.0477465im, 0.031831-0.063662im, 0.0795775-0.0795775im]))Since the structure of these targets can get complicated, e.g. nested tuples), the library also provides a set of functions for creating and resizing the velocity variable for in-place computations. For example:julia> vels = allocate_velocity(sys)\n(Complex{Float64}[0.0+0.0im, 0.0+0.0im, 0.0+0.0im], (Complex{Float64}[0.0+0.0im, 0.0+0.0im, 0.0+0.0im], Complex{Float64}[0.0+0.0im, 0.0+0.0im, 0.0+0.0im]))\n\njulia> induce_velocity!(vels, sys, src)\n(Complex{Float64}[0.0795775+0.0795775im, 0.031831+0.063662im, 0.0159155+0.0477465im], (Complex{Float64}[0.0+0.0im, -0.159155+0.0im, -0.0795775+0.0im], Complex{Float64}[0.0159155-0.0477465im, 0.031831-0.063662im, 0.0795775-0.0795775im]))The remaining sections of this page list the documentation for all the relevant methods for computing velocities. More detailed examples that show these methods working together can be found in the getting started guide and the Jupyter notebooks."
},

{
    "location": "velocities.html#VortexModel.Vortex.allocate_velocity",
    "page": "Computing Velocities",
    "title": "VortexModel.Vortex.allocate_velocity",
    "category": "Function",
    "text": "allocate_velocity(srcs)\n\nAllocate arrays of Complex128 to match the structure of srcs\n\nExample\n\njulia> points = Vortex.Point.(rand(Complex128, 2), rand(2));\n\njulia> blobs  = Vortex.Blob.(rand(Complex128, 3), rand(3), rand(3));\n\njulia> allocate_velocity(points)\n2-element Array{Complex{Float64},1}:\n 0.0+0.0im\n 0.0+0.0im\n\njulia> allocate_velocity((points, blobs))\n(Complex{Float64}[0.0+0.0im, 0.0+0.0im], Complex{Float64}[0.0+0.0im, 0.0+0.0im, 0.0+0.0im])\n\n\n\n"
},

{
    "location": "velocities.html#VortexModel.Vortex.reset_velocity!",
    "page": "Computing Velocities",
    "title": "VortexModel.Vortex.reset_velocity!",
    "category": "Function",
    "text": "reset_velocity!(vels[, srcs])\n\nSet all velocities in vels to zero\n\nIf srcs is provided, then the arrays in vels are resized their source counterpart, if necessary.\n\nExample\n\njulia> ẋs = (rand(Complex128, 1), rand(Complex128, 1))\n(Complex{Float64}[0.236033+0.346517im], Complex{Float64}[0.312707+0.00790928im])\n\njulia> points = Vortex.Point.(rand(Complex128, 2), rand(2));\n\njulia> blobs  = Vortex.Blob.(rand(Complex128, 3), rand(3), rand(3));\n\njulia> reset_velocity!(ẋs, (points, blobs));\n\njulia> ẋs\n(Complex{Float64}[0.0+0.0im, 0.0+0.0im], Complex{Float64}[0.0+0.0im, 0.0+0.0im, 0.0+0.0im])\n\n\n\n"
},

{
    "location": "velocities.html#VortexModel.Vortex.induce_velocity",
    "page": "Computing Velocities",
    "title": "VortexModel.Vortex.induce_velocity",
    "category": "Function",
    "text": "induce_velocity(target, element)\n\nCompute the velocity induced by element on target\n\ntarget can be:\n\na Complex128\na subtype of Vortex.PointSource\nan array or tuple of vortex elements\n\nwhile the element can be:\n\nany subtype of Vortex.Element\nan array or tuple of vortex elements\n\nExample\n\njulia> z = rand(Complex128)\n0.23603334566204692 + 0.34651701419196046im\n\njulia> point = Vortex.Point(z, rand());\n\njulia> srcs = Vortex.Point.(rand(Complex128, 10), rand(10));\n\njulia> induce_velocity(z, srcs[1])\n0.08722212007570912 + 0.14002850279102955im\n\njulia> induce_velocity(point, srcs[1])\n0.08722212007570912 + 0.14002850279102955im\n\njulia> induce_velocity(z, srcs)\n-0.4453372874427177 - 0.10592646656959151im\n\njulia> induce_velocity(point, srcs)\n-0.4453372874427177 - 0.10592646656959151im\n\n\n\n"
},

{
    "location": "velocities.html#VortexModel.Vortex.induce_velocity!",
    "page": "Computing Velocities",
    "title": "VortexModel.Vortex.induce_velocity!",
    "category": "Function",
    "text": "induce_velocity!(vels, target, element)\n\nCompute the velocity induced by element on target and store the result in vels\n\nvels should be the output of a call to allocate_velocity, target can be an array or tuple of vortex elements, while the element can be:\n\nany subtype of Vortex.Element\nan array or tuple of vortex elements\n\nExample\n\njulia> cluster₁ = Vortex.Point.(rand(Complex128, 5), rand(5));\n\njulia> cluster₂ = Vortex.Point.(rand(Complex128, 5), rand(5));\n\njulia> targets = (cluster₁, cluster₂);\n\njulia> sources = Vortex.Blob.(rand(Complex128), rand(10), 0.1);\n\njulia> ẋs = allocate_velocity(targets);\n\njulia> induce_velocity!(ẋs, targets, sources);\n\njulia> ẋs\n(Complex{Float64}[-1.28772-1.82158im, 1.9386-1.64147im, -1.56438+1.57158im, -0.626254+0.375842im, -0.806568-0.213201im], Complex{Float64}[-0.583672-2.26031im, -0.329778-1.43388im, 0.426927+1.55352im, -0.93755+0.241361im, -1.08949-0.35598im])\n\n\n\n"
},

{
    "location": "velocities.html#VortexModel.Vortex.self_induce_velocity!",
    "page": "Computing Velocities",
    "title": "VortexModel.Vortex.self_induce_velocity!",
    "category": "Function",
    "text": "self_induce_velocity!(vels, elements)\n\nCompute the self induced velocity of one or more vortex elements\n\nThis involves a recursive call to self_induce_velocity! and pairwise calls to mutually_induce_velocity!.\n\nExample\n\njulia> points = Vortex.Point.([-1, 1], 1.0)\n2-element Array{VortexModel.Vortex.Points.Point,1}:\n Point Vortex: z = -1.0 + 0.0im, Γ = 1.0\n Point Vortex: z = 1.0 + 0.0im, Γ = 1.0\n\njulia> vels = allocate_velocity(points)\n2-element Array{Complex{Float64},1}:\n 0.0+0.0im\n 0.0+0.0im\n\njulia> self_induce_velocity!(vels, points)\n\njulia> vels # should be ±0.25im/π\n2-element Array{Complex{Float64},1}:\n 0.0-0.0795775im\n 0.0+0.0795775im\n\n\n\n"
},

{
    "location": "velocities.html#VortexModel.Vortex.mutually_induce_velocity!",
    "page": "Computing Velocities",
    "title": "VortexModel.Vortex.mutually_induce_velocity!",
    "category": "Function",
    "text": "mutually_induce_velocity!(vs₁, vs₂, e₁, e₂)\n\nCompute the mutually induced velocities between e₁ and e₂ and store the results in vs₁ and vs₂\n\nThe default implementation simply calls induce_velocity! twice. This method is meant to be overwritten to take advantage of symmetries in certain pairwise vortex interations. For example, the velocity kernel for a point vortex is antisymmetric, so in computing the mutually induced velocities of two arrays of point vortices, we can half the number of calls to the velocity kernel.\n\n\n\n"
},

{
    "location": "velocities.html#Methods-1",
    "page": "Computing Velocities",
    "title": "Methods",
    "category": "section",
    "text": "allocate_velocity\nreset_velocity!\ninduce_velocity\ninduce_velocity!\nself_induce_velocity!\nmutually_induce_velocity!"
},

{
    "location": "velocities.html#Index-1",
    "page": "Computing Velocities",
    "title": "Index",
    "category": "section",
    "text": "Pages = [\"velocities.md\"]"
},

{
    "location": "timemarching.html#",
    "page": "Time Marching",
    "title": "Time Marching",
    "category": "page",
    "text": ""
},

{
    "location": "timemarching.html#Time-Marching-1",
    "page": "Time Marching",
    "title": "Time Marching",
    "category": "section",
    "text": "Coming soon..."
},

{
    "location": "noflowthrough.html#",
    "page": "Enforcing No-Flow-Through",
    "title": "Enforcing No-Flow-Through",
    "category": "page",
    "text": ""
},

{
    "location": "noflowthrough.html#Enforcing-No-Flow-Through-1",
    "page": "Enforcing No-Flow-Through",
    "title": "Enforcing No-Flow-Through",
    "category": "section",
    "text": "warning: Warning\nUnder construction...defddt1fracmathrmd1mathrmdt\n\nrenewcommandvecboldsymbol\nnewcommanduvec1vechat1\nnewcommandutangentuvectau\nnewcommandunormaluvecn\n\nrenewcommanddmathrmd\n\nnewcommandcrosstimes\nnewcommandabs1left1right\nnewcommandimmathrmi\nnewcommandeumathrme\nnewcommandpintint\nnewcommandconj11^star\nnewcommandRes2mathrmResleft(12right)\nnewcommandreal1mathrmReleft1right\nnewcommandimag1mathrmImleft1rightWe are interested in enforcing the no-flow-through condition on an infinitely thin, flat plate undergoing rigid body motion. The plate can be parameterized by its length, L, centroid position, vecc, and its angle of attack, alpha. Its motion is then specified by its centroid velocity, dotvecc, and angular velocity, dotvecalpha."
},

{
    "location": "noflowthrough.html#Vortex-Sheet-Strength-1",
    "page": "Enforcing No-Flow-Through",
    "title": "Vortex Sheet Strength",
    "category": "section",
    "text": "The plate is represented with a bound vortex sheet that constantly adjusts its circulation to enforce no-flow-through on its surface. We can show that the distribution of circulation, gamma, is governed by the following integral equation:<details>\n<summary></summary>\nThe no-flow-through condition requires that the component of fluid velocity normal to the sheet must be equal to the normal velocity of the sheet itself, i.e.\n$$\n\\begin{align*}\n    \\unormal \\cdot \\vec{u}(\\vec{x}_s)\n& = \\unormal \\cdot \\left[ \\dot{\\vec{c}} + \\dot{\\alpha} \\cross (\\vec{x}_s - \\vec{c}) \\right] \\\\\n& = \\left(\\unormal \\cdot \\vec{c}\\right) + \\dot{\\alpha} l\n\\end{align*}\n$$\nwhere\n\n- $\\vec{u}$ is the fluid velocity\n- $\\vec{x}_s$ is a position on the plate\n- $\\unormal$ is a unit vector normal to the plate\n- $l \\in [ -L/2, L/2 ] $ is distance between $\\vec{x}_s$ from the plate centroid\n\nWe can decompose the velocity field at any point in the fluid into contributions from the bound vortex sheet, $\\vec{u}_s$, and the free vorticity in the ambient fluid, $\\vec{u}_A$:\n$$\n\\vec{u}(\\vec{x}) = \\vec{u}_s(\\vec{x}) + \\vec{u}_A(\\vec{x}),\n$$\nso the no-flow-through condition can be written as:\n$$\n\\unormal \\cdot \\vec{u}_s(\\vec{x}) = \\left(\\unormal \\cdot \\vec{c}\\right) + \\dot{\\alpha} l - \\unormal \\cdot \\vec{u}_A(\\vec{x}).\n$$\n\nThe velocity field induced by a vortex sheet, $\\vec{u}_x(\\vec{x})$, is given by\n$$\n\\vec{u}_s(\\vec{x}) = \\frac{1}{2\\pi}\n\\int_\\mathcal{C} \\gamma(l) \\,\\uvec{k} \\cross\n\\frac{\\vec{x} - \\vec{x}_s(l)}{\\abs{\\vec{x} - \\vec{x}_s(l)}^2}\n\\d{l}\n$$\nwhere\n\n- $\\gamma$ is the strength of the sheet\n- $\\mathcal{C}$ is the curve occupied by the sheet\n- $\\uvec{k}$ is the unit vector point out of the plane.\n\nThe position along the vortex sheet can be expressed as\n$$\n\\vec{x}_s(l) = \\vec{c} + l\\utangent\n$$\nwhere $\\utangent$ is the unit tangent along the sheet.\nSimilarly, since we are interested in evaluating the velocity along the sheet, we can write\n$$\n\\vec{x}(l) = \\vec{c} + \\lambda\\utangent.\n$$\nWe can then write self-induced velocity of the bound vortex sheet as\n$$\n\\vec{u}_s(\\lambda) = \\frac{\\unormal}{2\\pi}\n\\int_{-\\frac{L}{2}}^\\frac{L}{2} \\frac{\\gamma(l)}{\\lambda - l}\n\\d{l}.\n$$\nSubstituting this expression back into the no-flow-through condition, we get\n</p>\n</details>beginequation\nfrac12pi\nint_-L2^L2 fracgamma(lambda)l - lambda\ndlambda\n= unormal cdot vecdotc\n+ dotalpha l\n- unormal cdot vecu_A(l)\nlabeleqintegral-equation\nendequationThe solution to this integral equation can be found in [Muskhelishvili]. If the velocity induced by ambient vorticity on the plate can be expanded into a Chebyshev series:unormal cdot vecu_Al(s) = sum_n = 0 A_n T_n(s)and Gamma_A is the total circulation in the ambient fluid, then the solution to eqrefeqintegral-equation can be written as:<details>\n<summary> </summary>\nTo make it easier to work with Chebyshev series, we will apply a change of variables $s := \\frac{2l}{L}$ so that the integral above goes from $-1$ to $1$:\n$$\n\\frac{1}{2\\pi}\n\\int_{-1}^1 \\frac{\\gamma(s)}{\\sigma - s}\n\\d{s}\n= \\unormal \\cdot \\vec{\\dot{c}}\n+ \\frac{\\dot{\\alpha}L}{2} \\sigma\n- \\unormal \\cdot \\vec{u}_A(\\sigma)\n$$\nFrom <a href=\"#footnote-Muskhelishvili\">[Muskhelishvili]</a>, we have that if\n$$\n\\frac{1}{\\pi\\im} \\int \\frac{\\varphi(t)}{t - t_0} \\d{t} = f(t_0)\n$$\nthen\n$$\n\\varphi(t_0) = \\frac{1}{\\pi\\im\\sqrt{t_0 - 1}\\sqrt{t_0 + 1}}\n\\int \\frac{\\sqrt{t - 1}\\sqrt{t + 1}}{t - t_0} f(t) \\d{t}\n+\n\\frac{P(t_0)}{\\sqrt{t_0 - 1}\\sqrt{t_0 + 1}}\n$$\nwhere $P$ is an arbitrary polynomial that must be chosen to satisfy far-field boundary conditions.\n\nIn our case, we have $\\varphi := \\im \\gamma$ and\n$$\nf := 2\\sum_{n = 0}^\\infty A_n T_n(\\sigma) - 2\\unormal \\cdot \\vec{\\dot{c}} - \\dot{\\alpha}L \\sigma\n$$\nso\n$$\n\\gamma(\\sigma)\n=\n\\frac{-2}{\\pi\\sqrt{1 - \\sigma}\\sqrt{1 + \\sigma}}\n\\int_{-1}^1 \\frac{\\sqrt{1 - s}\\sqrt{1 + s}}{s - \\sigma}\n\\left(\n\\sum_{n = 0}^\\infty A_n T_n(s) - \\unormal \\cdot \\vec{\\dot{c}} - \\frac{\\dot{\\alpha}L}{2} s\n\\right) \\d{s}\n+\n\\frac{P(t_0)}{\\sqrt{1 - \\sigma}\\sqrt{1 + \\sigma}}\n$$\n\nThe integral above is made of terms with the form\n$$\n\\pint_{-1}^1\n\\frac{\\sqrt{1 - s}\\sqrt{1 + s}}{s - \\sigma} T_n(s)\n\\d{s}\n$$\nwhich we can simplify using the properties of Chebyshev polynomials into\n$$\n\\pint_{-1}^1\n\\frac{\\sqrt{1 - s}\\sqrt{1 + s}}{s - \\sigma} T_n(s)\n\\d{s}\n=\n\\begin{cases}\n-\\pi T_1(\\sigma) & n = 0 \\\\\n-\\frac{\\pi}{2} T_2(\\sigma) & n = 1 \\\\\n-\\frac{\\pi}{2} \\left[T_{n+1}(\\sigma) - T_{n-1}(\\sigma)\\right] & n \\ge 2\n\\end{cases}.\n$$\nThis gives us\n$$\n\\gamma(\\sigma)\n=\n\\frac{-2}{\\pi\\sqrt{1 - \\sigma}\\sqrt{1 + \\sigma}}\n\\left\\{\n-\\pi A_0 \\sigma\n-\\frac{\\pi}{2} A_1\n+\\sum_{n = 1}^\\infty -\\frac{\\pi}{2}A_n \\left[T_{n+1}(\\sigma) - T_{n-1}(\\sigma)\\right]\n+ \\pi \\left(\\unormal \\cdot \\vec{\\dot{c}}\\right)\\sigma\n+ \\frac{\\pi}{2}T_2(\\sigma)\\frac{\\dot{\\alpha}L}{2}\n\\right\\}\n+\n\\frac{P(t_0)}{\\sqrt{1 - \\sigma}\\sqrt{1 + \\sigma}}.\n$$\n\nWe can find $P$ by satisfying Kelvin's circulation theorem.\nThis means that the amount of circulation contained in the bound vortex sheet should the negative of the circulation contained in the ambient vorticity, i.e.\n$$\n\\Gamma_s := \\int_{-\\frac{L}{2}}^{\\frac{L}{2}} \\gamma \\d{l} = -\\Gamma_A\n$$\n\nAgain, we use properties of Chebyshev polynomials to reduce the integral to\n$$\n\\begin{align*}\n\\frac{L}{2}\\int_{-1}^1 \\frac{P(s)}{\\sqrt{1 - s}\\sqrt{1 + s}} \\d{s} & = -\\Gamma_A,\n\\end{align*}\n$$\nwhich means that\n$$\nP = -\\frac{2\\Gamma_A}{L\\pi}.\n$$\n\nSo the final expression for the bound circulation is:\n</details>beginequation\ngammal(s) =\nfrac-frac2Gamma_ALpi + 2(A_0 - unormal cdot vecdotc) T_1(s) + (A_1 - fracdotalphaL2)T_2(s)sqrt1 - s^2 - 2sqrt1 - s^2sum_n = 2^infty A_n U_n-1(s)\nlabeleqgamma\nendequationnote: Note\nThis might look more similar to results from thin-airfoil theory if we rewrite the Chebyshev polynomials using trigonometric functions:gammal(theta) =\nfrac-frac2Gamma_ALpi + 2(A_0 - unormal cdot vecdotc) costheta + (A_1 - fracdotalphaL2)cos(2theta)sintheta - 2sum_n = 2^infty A_n sin(ntheta)The key difference is that we are free to relax the Kutta condition at the trailing edge."
},

{
    "location": "noflowthrough.html#Circulation-1",
    "page": "Enforcing No-Flow-Through",
    "title": "Circulation",
    "category": "section",
    "text": "In addition to the distribution of circulation along the plate, it will be useful to know the amount circulation contained between one end of the plate to an arbitrary point on its surface. By definition, we havebeginalign*\nGamma(l)  = int_-L2^l gamma(lambda) dlambda \nGammal(s) = fracL2int_-1^s gammal(sigma) dsigma\nendalign*We can integrate gamma term by term to obtain:<details>\n<summary></summary>\nIn equation $\\eqref{eq:gamma}$, the Chebyshev polynomial of the second kind in ther summation can be written in terms of Chebyshev polynomials of the first kind:\n$$\n2\\sqrt{1 - s^2}U_{n-1}(s)  = \\frac{T_{n-1}(s) - T_{n+1}(s)}{\\sqrt{1 - s^2}}.\n$$\nThis means that all the terms in equation $\\eqref{eq:gamma}$ can be expressed in the form:\n$$\n\\frac{T_n(s)}{\\sqrt{1 - s^2}}.\n$$\nThe integral of these terms are:\n$$\n\\begin{align*}\n\\int_{-1}^s \\frac{T_n(s)}{\\sqrt{1 - s^2}} \\d{s}\n& = \\int_{\\cos^{-1} s}^\\pi \\cos(n\\theta) \\d{\\theta} \\\\\n& = \\begin{cases}\n\\pi - \\cos^{-1}s &: n = 0 \\\\\n-\\frac{1}{n}\\sin\\left(n\\cos^{-1}s\\right) &: n > 0\n\\end{cases}.\n\\end{align*}\n$$\nWe can then multiply the expressions above with their corresponding coefficients to obtain:\n</details>Gammal(s)\n=Gamma_Aleft(fraccos^-1spi - 1right) - fracLsqrt1 - s^22left\n2left(A_0 - unormal cdot vecdotcright)\n+left(A_1 - fracdotalphaL2right)s\n+sum_n=2^infty\nA_n left(fracU_n(s)n+1 - fracU_n-2(s)n-1right)right[Muskhelishvili]: Muskhelishvili, Nikolaĭ Ivanovich, and Jens Rainer Maria Radok. Singular integral equations: boundary problems of function theory and their application to mathematical physics. Courier Corporation, 2008."
},

]}
