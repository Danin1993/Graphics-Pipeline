#include "driver_state.h"
#include <cstring>
#include <limits>
#include <vector>
#include <algorithm>

bool debug_mode = false;
bool its_clipping_time = true;

driver_state::driver_state()
{
}

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

    state.image_color = new pixel[width * height];
    state.image_depth = new float[width * height];
    for(size_t i = 0; i < width * height; i++) {
        state.image_depth[i] = std::numeric_limits<float>::max();
        state.image_color[i] = make_pixel(0, 0, 0);
    }
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
    const auto first_vertex = state.vertex_data;
    auto is_first_vertex = true;
    data_vertex in{};

    switch(type) {
        case render_type::triangle:
            for(int i = 0, j = 0; i < state.num_vertices; i++, j++) {
                tri[j].data = ptr;
                in.data = ptr;
                state.vertex_shader(in, tri[j], state.uniform_data);
                ptr += state.floats_per_vertex;
                if(j == 2) {
                    if(debug_mode)
                        std::cout << "Clipping triangle #" << i % 2 << std::endl;
                    if(its_clipping_time)
                        clip_triangle(state, (const data_geometry**) &tri, 0);
                    else
                        rasterize_triangle(state, (const data_geometry**) &tri);
                    j = -1;
                }
            }
            break;
        case render_type::indexed:
            for(int i = 0; i < 3 * state.num_triangles; i += 3) {
                for(int j = 0; j < 3; j++) {
                    ptr = &state.vertex_data[state.index_data[i + j] * state.floats_per_vertex];
                    tri[j].data = ptr;
                    in.data = ptr;
                    state.vertex_shader(in, tri[j], state.uniform_data);
                }
                if(debug_mode)
                    std::cout << "Clipping triangle #" << i % 2 << std::endl;
                if(its_clipping_time)
                    clip_triangle(state, (const data_geometry**) &tri, 0);
                else
                    rasterize_triangle(state, (const data_geometry**) &tri);
            }
            break;
        case render_type::fan:
            for(int i = 0, j = 0; i < state.num_vertices; i++, j++) {
                if(is_first_vertex) {
                    tri[j].data = first_vertex;
                    in.data = first_vertex;
                    state.vertex_shader(in, tri[j], state.uniform_data);
                    is_first_vertex = false;
                    i--;
                }
                else {
                    ptr += state.floats_per_vertex;
                    tri[j].data = ptr;
                    in.data = ptr;
                    state.vertex_shader(in, tri[j], state.uniform_data);
                }
                if(j == 2) {
                    if(debug_mode)
                        std::cout << "Clipping triangle #" << i % 2 << std::endl;
                    if(its_clipping_time)
                        clip_triangle(state, (const data_geometry**) &tri, 0);
                    else
                        rasterize_triangle(state, (const data_geometry**) &tri);
                    j = -1;
                    is_first_vertex = true;
                    ptr -= state.floats_per_vertex;
                    i--;
                }
            }
            break;
        case render_type::strip:
            for(int i = 0, j = 0; i < state.num_vertices; i++, j++) {
                tri[j].data = ptr;
                in.data = ptr;
                state.vertex_shader(in, tri[j], state.uniform_data);
                ptr += state.floats_per_vertex;
                if(j == 2) {
                    if(debug_mode)
                        std::cout << "Clipping triangle #" << i % 2 << std::endl;
                    if(its_clipping_time)
                        clip_triangle(state, (const data_geometry**) &tri, 0);
                    else
                        rasterize_triangle(state, (const data_geometry**) &tri);
                    j = -1;
                    ptr -= state.floats_per_vertex * 2;
                    i -= 2;
                }
            }
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
    if (face == 1) {
        rasterize_triangle(state, in);
        return;
    }
    else {
        vec4 A = (*in)[0].gl_Position;
        vec4 B = (*in)[1].gl_Position;
        vec4 C = (*in)[2].gl_Position;

        const data_geometry* new_in[3] = {in[0], in[1], in[2]};
        data_geometry first_dg[3];
        data_geometry second_dg[3];
        float new_A, new_B1, new_B2;
        vec4 P1, P2;

        // if all vertices of triangle are inside screen space, simply return
        if (A[2] < -A[3] && B[2] < -B[3] && C[2] < -C[3]) {
            return;
        } else {
            if (A[2] < -A[3] && B[2] >= -B[3] && C[2] >= -C[3]) {
                new_B1 = (-B[3] - B[2]) / (A[2] + A[3] - B[3] - B[2]);
                new_B2 = (-A[3] - A[2]) / (C[2] + C[3] - A[3] - A[2]);

                P1 = new_B1 * A + (1 - new_B1) * B;
                P2 = new_B2 * C + (1 - new_B2) * A;

                // =====================================================

                first_dg[0].data = new float[state.floats_per_vertex];
                first_dg[1] = *in[1];
                first_dg[2] = *in[2];

                for(int i = 0; i < state.floats_per_vertex; i++) {
                    switch (state.interp_rules[i]) {
                        case interp_type::flat:
                            first_dg[0].data[i] = (*in)[0].data[i];
                            break;
                        case interp_type::smooth:
                            first_dg[0].data[i] = new_B2 * (*in)[2].data[i] + (1 - new_B2) * (*in)[0].data[i];
                            break;
                        case interp_type::noperspective:
                            new_A = new_B2 * (*in)[2].gl_Position[3] / (new_B2 * (*in)[2].gl_Position[3] + (1 - new_B2) * (*in)[0].gl_Position[3]);
                            first_dg[0].data[i] = new_A * (*in)[2].data[i] + (1 - new_A) * (*in)[0].data[i];
                            break;
                        default:
                            break;
                    }
                }

                first_dg[0].gl_Position = P2;

                new_in[0] = &first_dg[0];
                new_in[1] = &first_dg[1];
                new_in[2] = &first_dg[2];

                clip_triangle(state, (const data_geometry**) &new_in, face + 1);

                // =====================================================

                second_dg[0].data = new float[state.floats_per_vertex];
                second_dg[1] = *in[1];
                second_dg[2] = *in[2];

                for (int i = 0; i < state.floats_per_vertex; ++i) {
                    switch (state.interp_rules[i]) {
                        case interp_type::flat:
                            second_dg[0].data[i] = (*in)[0].data[i];
                            break;
                        case interp_type::smooth:
                            second_dg[0].data[i] = new_B1 * (*in)[0].data[i] + (1 - new_B1) * (*in)[1].data[i];
                            break;
                        case interp_type::noperspective:
                            new_A = new_B1 * (*in)[0].gl_Position[3] / (new_B1 * (*in)[0].gl_Position[3] + (1 - new_B1) * (*in)[1].gl_Position[3]);
                            second_dg[0].data[i] = new_A * (*in)[0].data[i] + (1 - new_A) * (*in)[1].data[i];
                            break;
                        default:
                            break;
                    }
                }

                second_dg[0].gl_Position = P1;

                new_in[0] = &second_dg[0];
                new_in[1] = &first_dg[1];
                new_in[2] = &first_dg[0];
            }

            clip_triangle(state, (const data_geometry**) &new_in, face + 1);
        }
    }
}

// Rasterize the triangle defined by the three vertices in the "in" array.  This
// function is responsible for rasterization, interpolation of data to
// fragments, calling the fragment shader, and z-buffering.
void rasterize_triangle(driver_state& state, const data_geometry* in[3])
{
    int x[3], y[3], z[3];
    float w_div_2 = state.image_width / 2.0f;
    float h_div_2 = state.image_height / 2.0f;

    // calculates i and j coords in NDC for vertices
    for(int l = 0; l < 3; l++) {
        auto i = static_cast<int>(w_div_2 * ((*in)[l].gl_Position[0]/(*in)[l].gl_Position[3]) + (w_div_2 - 0.5f));
        auto j = static_cast<int>(h_div_2 * ((*in)[l].gl_Position[1]/(*in)[l].gl_Position[3]) + (h_div_2 - 0.5f));
        auto k = static_cast<int>(w_div_2 * ((*in)[l].gl_Position[2]/(*in)[l].gl_Position[3]) + (w_div_2 - 0.5f));
        x[l] = i;
        y[l] = j;
        z[l] = k;
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
    if(max_x > state.image_width - 1)
        max_x = state.image_width - 1;
    if(max_y > state.image_height - 1)
        max_y = state.image_height - 1;

    float area_ABC = (0.5f * ((x[1]*y[2] - x[2]*y[1]) - (x[0]*y[2] - x[2]*y[0]) + (x[0]*y[1] - x[1]*y[0])));

    auto *data = new float[MAX_FLOATS_PER_VERTEX];
    data_fragment frag_data{data};
    data_output output_data;

    // For each pixel in the bounding box of triangle, calculate it's barycentric weight with respect to the vertices
    // of the triangle. If pixel is in triangle, color it.
    for(int j = min_y + 1; j < max_y + 1; j++) {
        for(int i = min_x + 1; i < max_x + 1; i++) {
            float alpha_prime = (0.5f * ((x[1] * y[2] - x[2] * y[1]) + (y[1] - y[2])*i + (x[2] - x[1])*j)) / area_ABC;
            float beta_prime =  (0.5f * ((x[2] * y[0] - x[0] * y[2]) + (y[2] - y[0])*i + (x[0] - x[2])*j)) / area_ABC;
            float gamma_prime = (0.5f * ((x[0] * y[1] - x[1] * y[0]) + (y[0] - y[1])*i + (x[1] - x[0])*j)) / area_ABC;

            if (alpha_prime >= 0 && beta_prime >= 0 && gamma_prime >= 0) {
                float alpha = alpha_prime;
                float beta = beta_prime;
                float gamma = gamma_prime;

                float z_val = alpha * z[0] + beta * z[1] + gamma * z[2];

                // If z is closer than previously stored z, then color
                if(z_val < state.image_depth[i + j * state.image_width]) {
                    // update new z val
                    state.image_depth[i + j * state.image_width] = z_val;

                    for (int k = 0; k < state.floats_per_vertex; k++) {
                        float k_gour;
                        switch (state.interp_rules[k]) {
                            case interp_type::flat:
                                frag_data.data[k] = (*in)[0].data[k];
                                break;
                            case interp_type::smooth:
                                k_gour = (alpha / (*in)[0].gl_Position[3])
                                         + (beta / (*in)[1].gl_Position[3])
                                         + (gamma / (*in)[2].gl_Position[3]);

                                alpha_prime = alpha / (k_gour * (*in)[0].gl_Position[3]);
                                beta_prime = beta / (k_gour * (*in)[1].gl_Position[3]);
                                gamma_prime = gamma / (k_gour * (*in)[2].gl_Position[3]);
                            case interp_type::noperspective:
                                frag_data.data[k] = alpha_prime * (*in)[0].data[k]
                                                    + beta_prime * (*in)[1].data[k]
                                                    + gamma_prime * (*in)[2].data[k];
                                break;
                            default:
                                break;
                        }
                    }

                    state.fragment_shader(frag_data, output_data, state.uniform_data);

                    state.image_color[i + j * state.image_width] = make_pixel(
                            static_cast<int>(output_data.output_color[0] * 255),
                            static_cast<int>(output_data.output_color[1] * 255),
                            static_cast<int>(output_data.output_color[2] * 255));
                }
            }
        }
    }

    delete [] data;
}

data_geometry* create_triangle(driver_state& state, const data_geometry* in[3], vec4 A, vec4 B, vec4 C, int axis, int sign)
{
    float AB_t = ((sign * B[3] - B[axis]) / (A[axis] - sign*A[3] + sign*B[3] - B[axis]));
    float AC_t = ((sign * C[3] - C[axis]) / (A[axis] - sign*A[3] + sign*C[3] - C[axis]));
    float AB_k, AC_k;

    vec4 AB = AB_t * A + (1 - AB_t) * B;
    vec4 AC = AC_t * A + (1 - AC_t) * C;

    auto *new_tri = new data_geometry[3];
    new_tri[0].gl_Position = A;
    new_tri[1].gl_Position = AB;
    new_tri[2].gl_Position = AC;
    new_tri[0].data = new float[MAX_FLOATS_PER_VERTEX];
    new_tri[1].data = new float[MAX_FLOATS_PER_VERTEX];
    new_tri[2].data = new float[MAX_FLOATS_PER_VERTEX];

    for(int i = 0; i < state.floats_per_vertex; i++) {
        switch (state.interp_rules[i]) {
            case interp_type::flat:
                new_tri[0].data[i] = (*in)[0].data[i];
                break;
            case interp_type::smooth:
                new_tri[0].data[i] = (*in)[0].data[i];
                new_tri[1].data[i] = AB_t * (*in)[0].data[i] + (1-AB_t) * (*in)[1].data[i];
                new_tri[2].data[i] = AC_t * (*in)[0].data[i] + (1-AC_t) * (*in)[2].data[i];
                break;
            case interp_type::noperspective:
                AB_k = static_cast<float>(1.0 / (AB_t * A[3] + 1 - AB_t) * B[3]);
                AC_k = static_cast<float>(1.0 / (AC_t * A[3] + 1 - AC_t) * B[3]);
                AB_t *= A[3] * AB_k;
                AC_t *= A[3] * AC_k;
                new_tri[0].data[i] = (*in)[0].data[i];
                new_tri[1].data[i] = AB_t * (*in)[0].data[i] + (1-AB_t) * (*in)[1].data[i];
                new_tri[2].data[i] = AC_t * (*in)[0].data[i] + (1-AC_t) * (*in)[2].data[i];
                break;
            default:
                break;
        }
    }

    return new_tri;
}
