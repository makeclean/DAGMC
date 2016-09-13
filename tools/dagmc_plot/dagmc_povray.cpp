#include "moab/Core.hpp"
#include "moab/Types.hpp"
#include "MBTagConventions.hpp"
#include "moab/ProgOptions.hpp"

#include "uwuw.hpp"
#include "color_space.hpp"

#include <iostream>
#include <fstream>						
#include <set>

std::ofstream povout;

moab::Core *mbi;

/**
 * Generates a uniform color space for coloring material properties
 */
void generate_colors(int num_colors);

/**
 * Write the DAGMC geometry to a POV Mesh2 type input
 */
void write_mesh(std::vector<moab::EntityHandle> triangles, std::set<moab::EntityHandle> vertices,
		std::vector<double>,bool box = false);
/**
 * Adds a camera to the problem and therefore sets at which point we look at
 */
void add_camera(std::vector<double> origin, std::vector<double> look_at);

/**
 * Program to read a fully populated DAGMC Geometry and produce input for the 
 * POV-RAY rendering tool
 */
int main(int argc, char* argv[]) {

  moab::ErrorCode rval;

  std::string input_file;
  std::string output_file = "output.pov";

  std::vector<double> camera_origin; 
  std::vector<double> camera_look;
  std::vector<double> box_dims;

  bool box = false;

  // easy prog opts
  for ( int i = 1 ; i < argc ; i++ ) {
    if(std::string(argv[i]) == "-c"){
      camera_origin.push_back(std::atof(argv[i+1]));
      camera_origin.push_back(std::atof(argv[i+2]));
      camera_origin.push_back(std::atof(argv[i+3]));
    }

    if(std::string(argv[i]) == "-l"){
      camera_look.push_back(std::atof(argv[i+1]));
      camera_look.push_back(std::atof(argv[i+2]));
      camera_look.push_back(std::atof(argv[i+3]));
    }

    if(std::string(argv[i]) == "-b"){
      box_dims.push_back(std::atof(argv[i+1]));
      box_dims.push_back(std::atof(argv[i+2]));
      box_dims.push_back(std::atof(argv[i+3]));
      box_dims.push_back(std::atof(argv[i+4]));
      box_dims.push_back(std::atof(argv[i+5]));
      box_dims.push_back(std::atof(argv[i+6]));
    }

    if(std::string(argv[i]) == "-i"){
      input_file = std::string(argv[i+1]);
      std::ifstream f(input_file.c_str());
      if(!f.good()) {
	std::cout << "File doesnt exist" << std::endl;
	return -1;
      }
    }
    
    if(std::string(argv[i]) == "-o"){
      output_file = std::string(argv[i+1]);
    }
  }

  if (argc == 1 ) {
    std::cout << "Arguments needed!" << std::endl;
    return -1;
  }

  if(box_dims.size() == 6) {
    box = true;
  }

  if(camera_look.size() == 0 ) {
    std::cout << "No camera look set" << std::endl;
    return -1;
  }

  if(camera_origin.size() == 0 ) {
    std::cout << "No camera origin set" << std::endl;
    return -1;
  }

  if(input_file.length() == 0 ) {
    std::cout << "No filename set" << std::endl;
    return -1;
  }

  // open the output file
  povout.open(output_file);

  // try to load UWUW data to assign colours
  UWUW *uwuw = new UWUW(input_file);
  bool uwuw_data = (uwuw->material_library.size() > 0 ? true : false );
  if ( !uwuw_data ) {
    std::cout << "No materials founds, volumes will be coloured randomly" << std::endl;
  }

  colorspace *cs = new colorspace();
  std::vector<rgb_i> colorbar = cs->getRGBtable(uwuw->material_library.size());

  // new moab instance
  mbi = new moab::Core(); 
  // load a file
  rval = mbi->load_file(input_file.c_str());
  
  std::vector<moab::EntityHandle> meshsets;

  // the tags we need
  moab::Tag dim;
		     
  // need the dim tag in order to query the geometry type
  rval = mbi->tag_get_handle(GEOM_DIMENSION_TAG_NAME, 1, moab::MB_TYPE_INTEGER, dim, moab::MB_TAG_CREAT);

  // get all entity sets
  rval = mbi->get_entities_by_type(0,moab::MBENTITYSET,meshsets);
  int value = 0;

  // get the volumes
  std::vector<moab::EntityHandle> volumes;
  for ( unsigned int i = 0 ; i < meshsets.size() ; i++ ) {
    rval = mbi->tag_get_data(dim,&(meshsets[i]),1,&value);
    if ( value == 3 ) {
      volumes.push_back(meshsets[i]);
    }
  }

  // for each volume print out the surfaces 
  for ( unsigned int i = 0 ; i < volumes.size() ; i++ ) {
    std::vector<moab::EntityHandle> child_surfaces;
    rval = mbi->get_child_meshsets(volumes[i],child_surfaces);
    std::vector<moab::EntityHandle> triangles;
    // for each surface collect the triangles 
    for ( unsigned int j = 0 ; j < child_surfaces.size() ; j++ ) {
      std::vector<moab::EntityHandle> entities;
      rval = mbi->get_entities_by_type(child_surfaces[j],moab::MBTRI,entities);
      triangles.insert(triangles.begin(), entities.begin(), entities.end());
    }
    
    std::set<moab::EntityHandle> vertex_list;

    // for each triangle, collect the vertices into unique list
    for ( unsigned int j = 0 ; j < triangles.size() ; j++ ) {
      std::vector<moab::EntityHandle> vertices;
      rval = mbi->get_connectivity(&triangles[j],1,vertices);
      for ( unsigned int k = 0 ; k < vertices.size() ; k++ ) {
	       vertex_list.insert(vertices[k]);
      }
    }

    // now that we have the vertex data in a set we can determine 
    std::cout << "Writing volume " << i+1 << " of " << volumes.size() << std::endl;
    // get the color for this volume
    rgb_i color_value = colorbar[i];
    write_mesh(triangles,vertex_list,box_dims,box,color_value);
  }
  
  std::cout << "Adding camera ..." << std::endl;
  add_camera(camera_origin,camera_look);
  
  // close the output
  povout.close();
  return 0;
}

/**
 * Adds a camera to the POV-Ray file
 */
void add_camera(std::vector<double> origin, std::vector<double> look_at) {
  povout << "camera {" << std::endl;
  povout << "        location <";
  povout << origin[0] << ",";
  povout << origin[1] << ",";
  povout << origin[2] << ">" << std::endl;
  povout << "        look_at <";
  povout << look_at[0] << ",";
  povout << look_at[1] << ",";
  povout << look_at[2] << ">" << std::endl;
  povout << "}" << std::endl;
  return;
}

/**
 * writes a pov ray mesh object 
 */
void write_mesh(std::vector<moab::EntityHandle> triangles, std::set<moab::EntityHandle> vertices, 
		std::vector<double> box_dims, bool box, rgb_i color_value) {
  
  if ( box )
    povout << "difference {" << std::endl;

  povout << "mesh2 {" << std::endl;
  povout << "  vertex_vectors { " << std::endl;
  povout << "     " << vertices.size() << "," << std::endl;

  moab::ErrorCode rval;

  // make a local map of vertex id to 
  std::set<moab::EntityHandle>::iterator it;
  int id = 0;
  std::map<moab::EntityHandle,int> ent2id;
  for ( it = vertices.begin() ; it != vertices.end() ; ++it ) {
    // map of entity to id number 
    ent2id[*it] = id;
    //
    double coords[3];
    rval = mbi->get_coords(&(*it),1,coords);
    povout << "        ";
    povout << "<" << coords[0] << "," << coords[1];
    povout << "," << coords[2] << ">";
    if ( id == vertices.size()) 
      povout << std::endl;
    else
      povout << "," << std::endl;
    // povray meshes count from 0
    id++;
  }
  povout << "  }" << std::endl;

  // now make the triangles
  povout << "  face_indices {" << std::endl;
  povout << "   " << triangles.size() << "," << std::endl;
  for ( unsigned int i = 0 ; i < triangles.size() ; i++ ) {
    std::vector<moab::EntityHandle> vertices;
    rval = mbi->get_connectivity(&triangles[i],1,vertices);
    povout << "      <";
    povout << ent2id[vertices[0]] << ",";
    povout << ent2id[vertices[1]] << ",";
    povout << ent2id[vertices[2]] << ">";
    if( i != triangles.size() - 1 )
      povout << ",";
    povout << std::endl;
  }
  povout << "  }" << std::endl;

  // add a cutthrough box if desired
  if (box) {
    povout << "} " << std::endl;
    povout << "box { " << std::endl;
    povout << "<" << box_dims[0] << ",";
    povout << box_dims[1] << ",";
    povout << box_dims[2] << ">";
    povout << ",";
    povout << "<" << box_dims[3] << ",";
    povout << box_dims[4] << ",";
    povout << box_dims[5] << ">";
    povout << "} " << std::endl;
  }

  // give the true color
  povout << "pigment {color rgb<";
  povout << color_value.r << ",";
  povout << color_value.g << ",";
  povout << color_value.b << 
  povout << " >/255 }" << std::endl;
  povout << "}" << std::endl;
}
