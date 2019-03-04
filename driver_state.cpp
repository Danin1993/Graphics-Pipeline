#include "driver_state.h"
#include <cstring>

driver_state::driver_state() = default;

driver_state::~driver_state()
{
    delete [] image_color;
    delete [] image_depth;
}

// This function should allocate and initialize the arrays that store color and
// depth.  This is not done during the constructor since the width and height
// are not known when this class is constructed.
void initialize_render(driver_state& state, int width, int height)
{
    state.image_width=width;
    state.image_height=height;
    state.image_color=nullptr;
    state.image_depth=nullptr;
    //std::cout<<"TODO: allocate and initialize state.image_depth (when dealing with z-buffering)"<<std::endl;

    state.image_color = new pixel[width * height];
    for(size_t i = 0; i < width * height; i++)
        state.image_color[i] = make_pixel(0, 0, 0);
}

// This function will be called to render the data that has been stored in this class.
// Valid values of type are:
//   render_type::triangle - Each group of three vertices corresponds to a triangle.
//   render_type::indexed -  Each group of three indices in index_data corresponds
//                           to a triangle.  These numbers are indices into vertex_data.
//   render_type::fan -      The vertices are to be interpreted as a triangle fan.
//   render_type::strip -    The vertices are to be interpreted as a triangle strip.
void render(driver_state& state, render_type type)
{
    auto *tri = new data_geometry[3];
    auto ptr = state.vertex_data;
    data_vertex in{};

    switch(type) {
        case render_type::triangle:
            for(int i = 0, j = 0; i < state.num_vertices; i++, j++) {
                tri[j].data = ptr;
                in.data = ptr;
                state.vertex_shader(in, tri[j], state.uniform_data);
                ptr += state.floats_per_vertex;
                if(j == 2) {
                    rasterize_triangle(state, (const data_geometry**) &tri);
                    j = -1;
                }
            }
            break;
        case render_type::indexed:
            break;
        case render_type::fan:
            break;
        case render_type::strip:
            break;
        default:
            break;
    }

    delete [] tri;
}


// This function clips a triangle (defined by the three vertices in the "in" array).
// It will be called recursively, once for each clipping face (face=0, 1, ..., 5) to
// clip against each of the clipping faces in turn.  When face=6, clip_triangle should
// simply pass the call on to rasterize_triangle.
void clip_triangle(driver_state& state, const data_geometry* in[3],int face)
{
    if(face==6)
    {
        rasterize_triangle(state, in);
        return;
    }
    std::cout<<"TODO: implement clipping. (The current code passes the triangle through without clipping them.)"<<std::endl;
    clip_triangle(state,in,face+1);
}

// Rasterize the triangle defined by the three vertices in the "in" array.  This
// function is responsible for rasterization, interpolation of data to
// fragments, calling the fragment shader, and z-buffering.
void rasterize_triangle(driver_state& state, const data_geometry* in[3])
{
    int x[3], y[3];
    float w_div_2 = state.image_width / 2.0f;
    float h_div_2 = state.image_height / 2.0f;

    // calculates i and j coords in NDC for vertices
    for(int k = 0; k < 3; k++) {
        auto i = static_cast<int>(w_div_2 * ((*in)[k].gl_Position[0]/(*in)[k].gl_Position[3]) + (w_div_2 - 0.5f));
        auto j = static_cast<int>(h_div_2 * ((*in)[k].gl_Position[1]/(*in)[k].gl_Position[3]) + (h_div_2 - 0.5f));
        x[k] = i;
        y[k] = j;
    }

    // calculate min and max of triangle
    int min_x = std::min(std::min(x[0], x[1]), x[2]);
    int max_x = std::max(std::max(x[0], x[1]), x[2]);
    int min_y = std::min(std::min(y[0], y[1]), y[2]);
    int max_y = std::max(std::max(y[0], y[1]), y[2]);

    // check for cases where triangle goes off pixel grid
    if(min_x < 0)
        min_x = 0;
    if(min_y < 0)
        min_y = 0;
    if(max_x > state.image_width)
        max_x = state.image_width - 1;
    if(max_y > state.image_height)
        max_y = state.image_height - 1;

    float area_ABC = (0.5f * ((x[1]*y[2] - x[2]*y[1]) - (x[0]*y[2] - x[2]*y[0]) + (x[0]*y[1] - x[1]*y[0])));

    auto *data = new float[MAX_FLOATS_PER_VERTEX];
    data_fragment frag_data{data};
    //auto *frag_data = new data_fragment[MAX_FLOATS_PER_VERTEX];
    auto output_data = new data_output;

    // For each pixel in the bounding box of triangle, calculate it's barycentric weight with respect to the vertices
    // of the triangle. If pixel is in triangle, color it.
    for(int j = min_y; j < max_y + 1; j++) {
        for(int i = min_x; i < max_x + 1; i++) {
            float alpha = (0.5f * ((x[1] * y[2] - x[2] * y[1]) + (y[1] - y[2])*i + (x[2] - x[1])*j)) / area_ABC;
            float beta =  (0.5f * ((x[2] * y[0] - x[0] * y[2]) + (y[2] - y[0])*i + (x[0] - x[2])*j)) / area_ABC;
            float gamma = (0.5f * ((x[0] * y[1] - x[1] * y[0]) + (y[0] - y[1])*i + (x[1] - x[0])*j)) / area_ABC;

            if (alpha >= 0 && beta >= 0 && gamma >= 0) {
                for(int k = 0; k < state.floats_per_vertex; k++) {
                    switch(state.interp_rules[k]) {
                        case interp_type::flat:
                            frag_data.data[k] = in[0]->data[k];
                            break;
                        case interp_type::smooth:
                            break;
                        case interp_type::noperspective:
                            break;
                        default:
                            break;
                    }
                }

                state.fragment_shader(frag_data, *output_data, state.uniform_data);

                state.image_color[i + j * state.image_width] = make_pixel(
                        static_cast<int>(output_data->output_color[0] * 255),
                        static_cast<int>(output_data->output_color[1] * 255),
                        static_cast<int>(output_data->output_color[2] * 255));
            }
        }
    }

    delete [] data;
    delete output_data;
}