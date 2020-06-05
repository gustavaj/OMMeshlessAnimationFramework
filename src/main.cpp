/*
	Framework from here: https://github.com/SaschaWillems/Vulkan
	Some code for elephant including some shaders from here: https://prideout.net/blog/old/blog/index.html@tag=tessellation.html
	Some code from prototype
*/

#include "LatticeExample.cpp"

// TODO: Right now resources are not properly deleted when lattices are added from the menu

// TODO: Adding new local surface types are a lot of repeat work right now
// All the different pre-eval schemes are added in almost the same way but has to be done individually
// Maybe add some evaluator function that are used by all of them?

LatticeExample* example;

LRESULT CALLBACK WndProc(HWND hWnd, UINT uMsg, WPARAM wParam, LPARAM lParam)
{
	if (example != NULL)
	{
		example->handleMessages(hWnd, uMsg, wParam, lParam);
	}
	return (DefWindowProc(hWnd, uMsg, wParam, lParam));
}

int APIENTRY WinMain(HINSTANCE hInstance, HINSTANCE hPrevInstance, LPSTR pCmdLine, int nCmdShow)
{
	std::cerr << "asd" << std::endl;
	for (size_t i = 0; i < __argc; i++) { LatticeExample::args.push_back(__argv[i]); };

	bool debug = true;
	OML::LocalSurfaceType lsType = OML::LocalSurfaceType::Cubic_Bezier;
	OML::EvaluationMethod evalMethod = OML::EvaluationMethod::Direct;

	// Base example
	//example = new LatticeExample(debug, lsType, evalMethod);

	// Grid example
	example = new GridLatticeExample(debug, 500.0f, 500.0f, 10, 10, lsType, evalMethod);

	// Random grid example
	//example = new RandomGridLatticeExample(debug, 100.0f, 100.0f, 10, 10, lsType, evalMethod);

	// Cylinder Example
	//example = new CylinderLatticeExample(debug, 20.0f, 20.0f, 1, 4, lsType, evalMethod);

	// Sphere Example
	//example = new SphereLatticeExample(debug, 20.0f, 8, 8, lsType, evalMethod);

	// Torus Example
	//example = new TorusLatticeExample(debug, 500.0f, 250.0f, 6, 6, lsType, evalMethod);

	// Non-uniform grid example
	//example = new NonUniformGridExample(debug, lsType, evalMethod);

	// Grid where quads have angles != 90
	//example = new NonRectangularExample(debug, lsType, evalMethod);

	// T-locus example
	//example = new TLocusExample(debug, lsType, evalMethod);
	//example = new TLocusx4Example(debug, lsType, evalMethod);

	// Example with several lattices
	//example = new MultiLatticeExample(debug, lsType, evalMethod);

	// Example with all the evaluation methods for the same surface
	//example = new AllMethodsGridExample(debug, true, OML::LocalSurfaceType::Cubic_Bezier);

	// Example to test pixel acc in report.
	//example = new PixelAccTest(false);

	// Benchmarks
	//OML::LocalSurfaceType lsType = OML::LocalSurfaceType::Cubic_Bezier;
	//OML::EvaluationMethod evalMethod = OML::EvaluationMethod::Direct;
	//int rowcol = 60; // 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120
	//int tess = 16; // 16, 32, 64
	//example = new BenchmarkGrid(rowcol, tess, "./benchmarks/", lsType, evalMethod);

	OML::Timer::Start("init", "Setup example");

	example->initVulkan();
	example->setupWindow(hInstance, WndProc);
	example->prepare();

	OML::Timer::Stop("init", "Setup example");

	//example->lattices[0].setSurfaceColorMethod(OML::SurfaceColor::PixelAccuracy);
	//example->rebuildCmdBuffers = true;


	/*std::string base = "./benchmarks/direct_eval_bez3_500MB/";
	std::string filename = base + "tess16/" + std::to_string(rowcol) + "rc_lsType"
		+ std::to_string(static_cast<int>(lsType)) + "_eval" + std::to_string(static_cast<int>(evalMethod))
		+ "_2w_5d.txt";
	example->lattices[0].setTessellationFactors(16, 16);
	example->lattices[0].updateUniformBuffer();
	example->rebuildCmdBuffers = true;
	example->benchmark.filename = filename;
	example->renderLoop();


	filename = base + "tess32/" + std::to_string(rowcol) + "rc_lsType"
		+ std::to_string(static_cast<int>(lsType)) + "_eval" + std::to_string(static_cast<int>(evalMethod))
		+ "_2w_5d.txt";
	example->lattices[0].setTessellationFactors(32, 32);
	example->lattices[0].updateUniformBuffer();
	example->rebuildCmdBuffers = true;
	example->benchmark.filename = filename;
	example->renderLoop();

	filename = base + "tess64/" + std::to_string(rowcol) + "rc_lsType"
		+ std::to_string(static_cast<int>(lsType)) + "_eval" + std::to_string(static_cast<int>(evalMethod))
		+ "_2w_5d.txt";
	example->lattices[0].setTessellationFactors(64, 64);
	example->lattices[0].updateUniformBuffer();
	example->rebuildCmdBuffers = true;
	example->benchmark.filename = filename;*/
	example->renderLoop();




	delete(example);
}