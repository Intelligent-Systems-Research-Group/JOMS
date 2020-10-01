
#include "graph_def.h"
#include "MeshTopology.h"
#include "corpora.h"
#include "exporter.h"
#include <cxxopts.hpp>

int main(int argc, char* argv[]){   
    //int mode = 3;
    //writePersonFbx("", &mode, NULL, 0, 0, 0);
    
    cxxopts::Options options("trainer", "Model Trainer");
    options.add_options()
  	("o,out_path", "Output Directory", cxxopts::value<std::string>())
  	("s,script_path", "Opt Script Path", cxxopts::value<std::string>())
	("c,components", "Number of shape blend-shapes", cxxopts::value<int>())
	("d,dataset","Dataset ini Path", cxxopts::value<std::string>())
    ("m,model_weigths_path","Path to Model Weights HDF5 file", 
		cxxopts::value<std::string>()->default_value(""))
	("f,freeze_weights","Freeze Model Weights",
		cxxopts::value<bool>()->default_value("false")->implicit_value("true"))
	("n,npoints","Number of point per point cloud",cxxopts::value<int>())
	("p,dump_progress","Dump progress in out_path every p-th iteration",
		cxxopts::value<int>()->default_value("50"))
	("r,raycasting", "Iteration after which replace raycasting with icp",
		cxxopts::value<int>()->default_value("180"))
	("j,jointly","Iteration after which we jointly optimize in R^n and surface correspondences",
		cxxopts::value<int>()->default_value("180"))
	("a,alternating",
		"Iteration interval after which to interleave discrete correspondences optimization",
		cxxopts::value<int>()->default_value("5"))
	("b,base_dataset_folder","path to dataset folder",cxxopts::value<std::string>())
	("u,use_model_to_data","use model to data error term",
		cxxopts::value<bool>()->default_value("false")->implicit_value("true"))
	("D,deform_joints","number of muscle deform joints influence",cxxopts::value<int>())
	("i,iteration_data","Iteration after which the point cloud data is used",cxxopts::value<int>())
	("w,continuous_step_size","Continuous step size",cxxopts::value<float>())
	("e,iteration_end","Number of outer LM Iterations",cxxopts::value<int>()->default_value("999"))
	("t,temporal_regularizer","Smooth pose regularizer",
   		 cxxopts::value<bool>()->default_value("false")->implicit_value("true"))
	("T,template_path","path to template folder",cxxopts::value<std::string>()->default_value("template"))
	("J,njoints", "Number of joints of the skeleton", cxxopts::value<int>())
	("V,max_vertex_ring","Max vertices on joint boundaries", cxxopts::value<int>())
	("C,center_vertex_id","Center Vertex Id",cxxopts::value<int>())
	("S,enable_symmetry","Enable symmetry",cxxopts::value<bool>()->default_value("false")->implicit_value("true"))
	("h,help", "Print usage")
  ;
    auto result = options.parse(argc, argv);
	if (result.count("help"))
    {
      std::cout << options.help() << std::endl;
      exit(0);
    }
    std::string dpath = result["dataset"].as<std::string>();

    Corpora corpora = constructFromConfig(dpath);

    //HeapProfilerStart();
    InputParams params;
	params.center_vertex_id = result["center_vertex_id"].as<int>();
	params.template_path = result["template_path"].as<std::string>();
    params.template_data_marker_path = params.template_path + "/markers_data2.txt";
    params.template_marker_path = params.template_path + "/markersE2.txt";
    params.template_vgroups_path = params.template_path + "/vgroupsE2.txt";
    params.template_mesh_path = params.template_path + "/shape.obj";
    params.template_laplace_path = params.template_path + "/laplacian.txt";
    params.template_skeleton_path = params.template_path + "/skeleton.json";
	params.template_surface_map_path = params.template_path + "/surface_to_keypoints.json";

	params.use_symmetry = result["enable_symmetry"].as<bool>();
    std::map<std::string, std::vector<size_t>> marker_vertex_groups;
#ifdef USE_MARKERS
	//loadVertexGroups("markers_fixed.txt", &marker_vertex_groups, NULL, NULL); //"markers.txt.old"
    loadVertexGroups(params.template_marker_path, &marker_vertex_groups, NULL, NULL);
#endif
    //@TODO To config
	params.iteration_data = result["iteration_data"].as<int>();
	params.useModel2Data = result["use_model_to_data"].as<bool>();
	params.base_dataset_folder = result["base_dataset_folder"].as<std::string>();
	params.raycasting = result["raycasting"].as<int>();
	params.jointly = result["jointly"].as<int>();
	params.alternating = result["alternating"].as<int>();
	params.dump_progress = result["dump_progress"].as<int>();
	params.freeze_weights =  result["freeze_weights"].as<bool>();
	params.model_weights_path = result["model_weigths_path"].as<std::string>();
    params.optscript_path = result["script_path"].as<std::string>();
    params.out_path = result["out_path"].as<std::string>();
    params.ncomp = result["components"].as<int>();
    params.npoints = result["npoints"].as<int>(); //6000; //40000 6890 6000
    params.nlabels = marker_vertex_groups.size();//marker_vertex_groups.size(); //100
	params.max_vertex_ring = result["max_vertex_ring"].as<int>();
	{
		MeshTopology top;
		Matrix3X control_vertices;
	    loadObj(params.template_mesh_path, &top, &control_vertices,NULL,1.0f);    
   		params.nverts = control_vertices.cols();
    }
    params.njoints = result["njoints"].as<int>(); //33 17
	params.njointlabels = params.njoints;
    params.nvarjoints = params.njoints; //14 33
    params.ndeformjoints = result["deform_joints"].as<int>();
    params.ndeformshapes = 9;
    params.nbasis = 1;
	params.continuous_step_size = result["continuous_step_size"].as<float>();
	params.iteration_end = result["iteration_end"].as<int>();    
    params.nscans = corpora.scan_names.size();
    params.npersons = corpora.person_ids.size(); // 512

	params.camPaths = {params.template_path+"/cameras.json",params.template_path+"/caesar_cameras.json"};
	params.nviews = 8;
	params.useTemporal = result["temporal_regularizer"].as<bool>();
	//params.camPaths = {"template/cameras_faust.json"};
	//params.nviews = 10;

	std::cout << "Freeze is : " << params.freeze_weights << std::endl;
    Problem prob(params, corpora);
    prob.run();
    //HeapProfilerStop();
}
