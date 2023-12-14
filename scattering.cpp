#include <iostream>
#include <cmath>
#include <random>
#include <fstream>
#include <algorithm>
#include <iterator>
#include <stdlib.h>
#include <sys/stat.h>
#include <vector>
#include <iomanip>
#include <sstream>

using namespace std;

float const impulse_coefficient = 0.000126;
int const cell_size = 200;

struct Atom
{
    float* impulse;
    float* coordinate_0;
    float* coordinate;
    unsigned int photon_absorption;
    unsigned int collisions_laser;
    unsigned int collisions_chaos;
    unsigned int photons_inside;
    unsigned int photons_tmp;
    float radius;
    bool tmp;
    Atom(float* impulse, float* coordinate_0, float* coordinate, 
        unsigned int photon_absorption, unsigned int collisions_laser, unsigned int collisions_chaos, // photon absorption - change cell
        unsigned int photons_inside, unsigned int photons_tmp, float radius, bool tmp) :
            impulse{ impulse }, coordinate_0{ coordinate_0 }, coordinate{ coordinate }, photon_absorption{ photon_absorption }, collisions_laser{collisions_laser}, collisions_chaos{ collisions_chaos },
            photons_inside{ photons_inside }, photons_tmp{ photons_tmp }, radius{radius}, tmp{ tmp} { }
    Atom():
        impulse{new float[3] {0., 0., 0.}}, coordinate_0{ new float[3] {0., 0., 0.}}, coordinate{ new float[3] {0., 0., 0.}}, photon_absorption{0 }, collisions_laser{0}, collisions_chaos{ 0 },
        photons_inside{ 0 }, photons_tmp{ 0 }, radius{0.24}, tmp{ false} { }
    void print(){
        std::cout << coordinate[0] << ' ' << coordinate[1] << ' ' << coordinate[2] << ' ' << photons_inside << ' ' << photons_tmp << ' ' << tmp << endl;
    }
    void write(std::ofstream& out){
        out << coordinate_0[0] << ' ' << coordinate_0[1] << ' ' << coordinate_0[2] << ' ';
        out << impulse[0] << ' ' << impulse[1] << ' ' << impulse[2] << ' ';
        out << coordinate[0] << ' ' << coordinate[1] << ' ' << coordinate[2] << ' ';
        out << collisions_laser << ' ' << collisions_chaos << ' ' ;
        out << photons_inside << ' ' << photons_tmp << ' ' << tmp << std::endl;

    }
    void move(){
        coordinate[0] += impulse[0];
        coordinate[1] += impulse[1];
        coordinate[2] += impulse[2];
    }
};

struct Photon
{
    float* emission;
    float* dir;
    Photon(float* emission, float* dir): emission{emission}, dir{dir} {}
    void print(){
        std::cout << dir[0] << ' ' << dir[1] << ' ' << dir[2] << ' ' << emission[0] << ' ' << emission[1] << ' ' << emission[2] << endl;
    }
};

struct Ball_cell
{
    float* coordinate_centre;
    Atom* elements;
    unsigned int number_of_elements;
    float edge;
    Ball_cell():
        coordinate_centre{new float[3]}, elements{new Atom[cell_size]}, number_of_elements{ 0 }, edge{ 0 } {}
    void print(){ 
        std::cout << coordinate_centre[0] << ' ' << coordinate_centre[1] << ' ' << coordinate_centre[2] << ' ' << edge << ' ' << number_of_elements << endl;
    }
};

void write(std::string file, Ball_cell* ball_cells, int const ball_cells_size, Atom* out_of_cloud, int* size_out_of_cloud){
    std::ofstream out;
    out.open(file);
    if (out.is_open()) {
        for (int i = 0; i < ball_cells_size; i++){
            for (int j = 0; j < ball_cells[i].number_of_elements; j++){
                if (ball_cells[i].elements[j].tmp == false){
                    continue;
                }
                ball_cells[i].elements[j].write(out);
            }
        }

        for (int i = 0; i < *size_out_of_cloud; i++){
            if (out_of_cloud[i].tmp == false){
                cout << 539 << endl;
            }
            out_of_cloud[i].write(out);
        }
    }
    out.close();    
    return;
}

float distance_between_points(float* point_1, float* point_2){
    return sqrt((point_1[0] - point_2[0])*(point_1[0] - point_2[0]) + (point_1[1] - point_2[1])*(point_1[1] - point_2[1]) + (point_1[2] - point_2[2])*(point_1[2] - point_2[2]));
}

float distance_between_line_and_point(float* dir, float* p_1, float* p_2){
    return sqrt(((dir[1]*(p_1[2] - p_2[2])-dir[2]*(p_1[1]-p_2[1])) * (dir[1]*(p_1[2] - p_2[2])-dir[2]*(p_1[1]-p_2[1]))
                    + (dir[0]*(p_1[2] - p_2[2])-dir[2]*(p_1[0]-p_2[0])) * (dir[0]*(p_1[2] - p_2[2])-dir[2]*(p_1[0]-p_2[0])) 
                    + ((dir[0]*(p_1[1] - p_2[1])-dir[1]*(p_1[0]-p_2[0]))) * ((dir[0]*(p_1[1] - p_2[1])-dir[1]*(p_1[0]-p_2[0]))))/(dir[0]*dir[0] + dir[1]*dir[1] + dir[2]*dir[2]));
}

float scalar_product(float* a, float *b){
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

void collision(Photon* photon, Atom* atom){
    atom->photons_tmp += 1;
    if (photon->dir[0] == 1){ 
        atom->collisions_laser += 1;
    }
    else{
        atom->collisions_chaos += 1;
    }

    atom->impulse[0] += photon->dir[0] * impulse_coefficient;
    atom->impulse[1] += photon->dir[1] * impulse_coefficient;
    atom->impulse[2] += photon->dir[2] * impulse_coefficient;

    if (atom->photons_tmp > 2) std::cout << 98 << endl;
}

void emission(Photon* photon, Atom* atom){
    atom->photons_tmp -= 1;
    
    atom->impulse[0] -= photon->dir[0] * impulse_coefficient;
    atom->impulse[1] -= photon->dir[1] * impulse_coefficient;
    atom->impulse[2] -= photon->dir[2] * impulse_coefficient;
    
}

void fly_away(Photon* photon, unsigned int number_of_photons, unsigned long long* photon_out, int const screen_size, int const screen_resolution, vector< vector<int> >* screen, int const screen_len){
    (*photon_out) += number_of_photons;
    if (photon->dir[0] == 1. ){
        if (photon->emission[1] + screen_size < 0 || photon->emission[1] + screen_size > 2 * screen_size){
            photon->print();
            cout << 147 << endl;
            
        }
        if (photon->emission[2] + screen_size < 0 || photon->emission[2] + screen_size > 2 * screen_size){
            photon->print();
            cout << 153 << endl;
            
        }

        (*screen)[static_cast<int>((photon->emission[1] + screen_size) * screen_resolution)][static_cast<int>((photon->emission[2] + screen_size) * screen_resolution)] += number_of_photons;
    }

}

Atom* find_nearest_atom(Photon* photon, Atom* atom_departure, Ball_cell* ball_cells, int ball_cells_size){
    Atom* nearest_atom = nullptr;
    for (int i = 0; i < ball_cells_size; i++){
        if(ball_cells[i].number_of_elements == 0){
            continue;
        }

        if ((photon->dir[0]*(ball_cells[i].coordinate_centre[0] - photon->emission[0]) + //TODO разность массивов
            photon->dir[1]*(ball_cells[i].coordinate_centre[1] - photon->emission[1]) + 
            photon->dir[2]*(ball_cells[i].coordinate_centre[2] - photon->emission[2])) < 0
            and distance_between_points(atom_departure->coordinate, ball_cells[i].coordinate_centre) > 1.5*ball_cells[i].edge){
                continue;
        }
        if (nearest_atom != nullptr){
            if (distance_between_points(nearest_atom->coordinate, atom_departure->coordinate) + ball_cells[i].edge/2 <\
                distance_between_points(atom_departure->coordinate, ball_cells[i].coordinate_centre)){
                    continue;
            }
        }

        if (distance_between_line_and_point(photon->dir, photon->emission, ball_cells[i].coordinate_centre) < 1.72*ball_cells[i].edge/2){
            
            for (int j = 0; j < ball_cells[i].number_of_elements; j++){
                if (ball_cells[i].elements[j].tmp == false){
                    continue;
                }

                if (&ball_cells[i].elements[j] == atom_departure){
                    continue;
                }

                if (photon->dir[0]*(ball_cells[i].elements[j].coordinate[0] - photon->emission[0]) + 
                    photon->dir[1]*(ball_cells[i].elements[j].coordinate[1] - photon->emission[1]) + 
                    photon->dir[2]*(ball_cells[i].elements[j].coordinate[2] - photon->emission[2]) < 0){
                        continue;
                }

                if (distance_between_line_and_point(photon->dir, photon->emission, ball_cells[i].elements[j].coordinate) < ball_cells[i].elements[j].radius){  // 0.24 = (sigma_0/pi)**0.5/sigma_r*r
                    if (not nearest_atom){
                        nearest_atom = &ball_cells[i].elements[j];
                    }
                    else if (distance_between_points(ball_cells[i].elements[j].coordinate, atom_departure->coordinate) < distance_between_points(nearest_atom->coordinate, atom_departure->coordinate)){
                        nearest_atom = &ball_cells[i].elements[j];
                    }
                }
            }
        }
    }
    return nearest_atom;
}

void next_collision(Photon* photon, Atom* atom_departure, int number_of_photons, Ball_cell* ball_cells, int ball_cells_size, unsigned long long* photon_out, int const screen_size, int const screen_resolution, vector< vector<int> >* screen, int const screen_len){
    if (number_of_photons == 1){
        Atom* nearest_atom = find_nearest_atom(photon, atom_departure, ball_cells, ball_cells_size);
       
        if (nearest_atom){
            collision(photon, nearest_atom);
            if (nearest_atom->photons_tmp > 1){

                emission(photon, nearest_atom);
                emission(photon, nearest_atom);
                std::copy_n(nearest_atom->coordinate, 3, photon->emission);
                

                Atom* next_atom = find_nearest_atom(photon, nearest_atom, ball_cells, ball_cells_size);

                if (not next_atom){
                    fly_away(photon, 2, photon_out, screen_size, screen_resolution, screen, screen_len);
                    return;
                }

                else if (next_atom->photons_tmp == 0){
                    collision(photon, next_atom);
                    std::copy_n(next_atom->coordinate, 3, photon->emission);
                    next_collision(photon, next_atom, 1, ball_cells, ball_cells_size, photon_out, screen_size, screen_resolution, screen, screen_len);
                    return;
                }

                else if (next_atom->photons_tmp == 1){
                    std::copy_n(next_atom->coordinate, 3, photon->emission);
                    emission(photon, next_atom);
                    next_collision(photon, next_atom, 3, ball_cells, ball_cells_size, photon_out, screen_size, screen_resolution, screen, screen_len);
                    return;
                }
            }
            else {
                return;
            }
        }
        else{
            fly_away(photon, 1, photon_out, screen_size, screen_resolution, screen, screen_len);
            return;
        }
        return;
    }

    else if (number_of_photons >= 2){
        Atom* nearest_atom = find_nearest_atom(photon, atom_departure, ball_cells, ball_cells_size);
        if (nearest_atom != nullptr){
            std::copy_n(nearest_atom->coordinate, 3, photon->emission);
            if (nearest_atom->photons_tmp == 0){
                collision(photon, nearest_atom);
                next_collision(photon, nearest_atom, number_of_photons - 1, ball_cells, ball_cells_size, photon_out, screen_size, screen_resolution, screen, screen_len);
                return;
            }
            
            if (nearest_atom->photons_tmp == 1){
                
                emission(photon, nearest_atom);
                Atom* next_atom = find_nearest_atom(photon, nearest_atom, ball_cells, ball_cells_size);

                if (not next_atom){
                    fly_away(photon, number_of_photons+1, photon_out,screen_size, screen_resolution, screen, screen_len);
                    return;
                }
                
                else if (next_atom->photons_tmp == 0){
                    collision(photon, next_atom);
                    std::copy_n(next_atom->coordinate, 3, photon->emission);
                    next_collision(photon, next_atom, number_of_photons, ball_cells, ball_cells_size, photon_out, screen_size, screen_resolution, screen, screen_len);
                    return;
                }

                else if (next_atom->photons_tmp == 1){
                    emission(photon, next_atom);
                    std::copy_n(next_atom->coordinate, 3, photon->emission);
                    next_collision(photon, next_atom, number_of_photons+2, ball_cells, ball_cells_size, photon_out, screen_size, screen_resolution, screen, screen_len);
                    return;
                }
            }
        }
        else{
            fly_away(photon, number_of_photons, photon_out, screen_size, screen_resolution, screen, screen_len);
            return;
        }
        return;
    }
}

void emit_photon(Photon* photon, Atom* atom, float* random_array, int random_array_size, unsigned long long* counter){
    

    float phi = random_array[*counter]*2*3.14;
    (*counter)++;

    if (*counter >= random_array_size){
        *counter = rand()%(random_array_size-5);
    }
    
    float theta = random_array[*counter]*2*3.14;
    (*counter)++;

    if (*counter >= random_array_size){
        *counter = rand()%(random_array_size-5);
    }   
    
    std::copy_n(atom->coordinate, 3, photon->emission);
    
    photon->dir[0] = cos(theta) * cos(phi);
    photon->dir[1] = cos(theta) * sin(phi);
    photon->dir[2] = sin(theta);
    
    atom->impulse[0] -= photon->dir[0] * impulse_coefficient;
    atom->impulse[1] -= photon->dir[1] * impulse_coefficient;
    atom->impulse[2] -= photon->dir[2] * impulse_coefficient;
}

void change_cell(Atom* atom, Ball_cell* cell_init, Ball_cell* ball_cells, int ball_cells_size, unsigned long long* photon_out, Atom* out_of_cloud, int* out_of_cloud_size, unsigned long long* photon_in){

    if (distance_between_points(atom->coordinate, cell_init->coordinate_centre) < sqrt(3)*cell_init->edge/2){
        return;
    }

    atom->tmp = false;
    
    bool break_flag = false;
    for (int i = 0; i < ball_cells_size; i++){
        if (distance_between_points(atom->coordinate, ball_cells[i].coordinate_centre) < sqrt(3)*ball_cells[i].edge/2){
            for (int j = 0; j < ball_cells[i].number_of_elements+2; j++){
                if (ball_cells[i].elements[j].tmp == false){
                    //ball_cells[i].elements[j] = *atom; // TODO Оператор присваивания

                    ball_cells[i].elements[j].tmp = true;
                    std::copy_n(atom->coordinate_0, 3, ball_cells[i].elements[j].coordinate_0);
                    std::copy_n(atom->coordinate, 3, ball_cells[i].elements[j].coordinate);
                    std::copy_n(atom->impulse, 3, ball_cells[i].elements[j].impulse);
                    ball_cells[i].elements[j].photons_inside = atom->photons_inside;
                    ball_cells[i].elements[j].photons_tmp = atom->photons_tmp;
                    ball_cells[i].elements[j].radius = atom->radius;
                    ball_cells[i].elements[j].photon_absorption = atom->photon_absorption;
                    ball_cells[i].elements[j].collisions_chaos = atom->collisions_chaos;
                    ball_cells[i].elements[j].collisions_laser = atom->collisions_laser;
                    
                    break_flag = true;
                    
                    if (j >= ball_cells[i].number_of_elements){
                        ball_cells[i].number_of_elements += 1;
                    }
                    break;
                }
            }
        }
        if (break_flag){
            break;
        }
    }

    if (break_flag == false){
        
        if (atom->photons_tmp == 1 || atom->photons_inside == 1){
            (*photon_out)++;
            atom->photons_tmp = 0;
            atom->photons_inside = 0;
        }
        out_of_cloud[*out_of_cloud_size].tmp = true;
        std::copy_n(atom->coordinate_0, 3, out_of_cloud[*out_of_cloud_size].coordinate_0);
        std::copy_n(atom->coordinate, 3, out_of_cloud[*out_of_cloud_size].coordinate);
        std::copy_n(atom->impulse, 3, out_of_cloud[*out_of_cloud_size].impulse);
        out_of_cloud[*out_of_cloud_size].photons_inside = atom->photons_inside;
        out_of_cloud[*out_of_cloud_size].photons_tmp = atom->photons_tmp;
        out_of_cloud[*out_of_cloud_size].radius = atom->radius;
        out_of_cloud[*out_of_cloud_size].photon_absorption = atom->photon_absorption;
        out_of_cloud[*out_of_cloud_size].collisions_chaos = atom->collisions_chaos;
        out_of_cloud[*out_of_cloud_size].collisions_laser = atom->collisions_laser;

        (*out_of_cloud_size)++;  

        if (*out_of_cloud_size >= 10000){
            cout << 360 << endl;
        }
    }    
}

void clear(Ball_cell* ball_cells, int const ball_cells_size, unsigned long long* photon_out, unsigned long long* photon_in, Atom* out_of_cloud, int* out_of_cloud_size){
    *photon_in = 0;
    *photon_out = 0;
    for (int i = 0; i < ball_cells_size; i++){
        for (int j = 0; j < ball_cells[i].number_of_elements; j++){
            if (ball_cells[i].elements[j].tmp == 0){
                continue;
            }
            std::copy_n(ball_cells[i].elements[j].coordinate_0, 3, ball_cells[i].elements[j].coordinate);
            ball_cells[i].elements[j].photons_inside = 0;
            ball_cells[i].elements[j].photons_tmp = 0;
            ball_cells[i].elements[j].impulse[0] = 0.;
            ball_cells[i].elements[j].impulse[1] = 0.;
            ball_cells[i].elements[j].impulse[2] = 0.;
            ball_cells[i].elements[j].collisions_chaos = 0;
            ball_cells[i].elements[j].collisions_laser = 0;
            ball_cells[i].elements[j].photon_absorption = 0;
            change_cell(&ball_cells[i].elements[j], &ball_cells[i], ball_cells, ball_cells_size, photon_out, out_of_cloud, out_of_cloud_size, photon_in);
        }
    }

    bool break_flag = false;
    for (int i = 0; i < *out_of_cloud_size; i++){
        out_of_cloud[i].tmp = 0;

        break_flag = false;
        for (int j = 0; j < ball_cells_size; j++){
            if (distance_between_points(out_of_cloud[i].coordinate, ball_cells[j].coordinate_centre) <  1.72 * ball_cells[j].edge/2){
                for(int k = 0; k < ball_cells[j].number_of_elements + 3; k++){
                    if (ball_cells[j].elements[k].tmp == 0){
                        std::copy_n(out_of_cloud[i].coordinate_0, 3, ball_cells[j].elements[k].coordinate);
                        std::copy_n(out_of_cloud[i].coordinate_0, 3, ball_cells[j].elements[k].coordinate_0);
                        ball_cells[j].elements[k].photons_inside = 0;
                        ball_cells[j].elements[k].impulse[0] = 0.;
                        ball_cells[j].elements[k].impulse[1] = 0.;
                        ball_cells[j].elements[k].impulse[2] = 0.;
                        ball_cells[j].elements[k].collisions_chaos = 0;
                        ball_cells[j].elements[k].collisions_laser = 0;
                        ball_cells[j].elements[k].photon_absorption = 0; 
                        ball_cells[j].elements[k].tmp = 1;
                        break_flag = true;
                        break;
                    }
                    
                }
                if (not break_flag){
                    cout << 450 << endl;
                }

                if (ball_cells[j].number_of_elements > cell_size - 2){
                    cout << 455 << endl;
                }
                

            }
            if (break_flag){
                break;
            }
        }
    }
    *out_of_cloud_size = 0;

    int number_of_atoms = 0;
    for (int i = 0; i < ball_cells_size; i++){
        for (int j = 0; j < ball_cells[i].number_of_elements; j++){
            if (ball_cells[i].elements[j].tmp == 1){
                number_of_atoms++;
            }
        }
    }
    std::cout << number_of_atoms << endl;

}

void clear_cells(Ball_cell* ball_cells, int const ball_cells_size){ //TODO написать функцию, которая будет очищать клетки
    for (int i = 0; i < ball_cells_size; i++){
        if (ball_cells[i].number_of_elements == 0){
            continue;
        }
        while (ball_cells[i].elements[ball_cells[i].number_of_elements-1].tmp != 1){
            ball_cells[i].number_of_elements --;
            if (ball_cells[i].number_of_elements == 0){
                break;
            }
        }
    }
}

void experiment(float intensity_coeff, Ball_cell* ball_cells, int const ball_cells_size, Atom* out_of_cloud, int* size_out_of_cloud){

    int const screen_size = 56;
    //Детектор фотонов
    int const screen_resolution = 5;
    int const screen_len = (screen_size*2*screen_resolution) * (screen_size*2*screen_resolution);
    vector <vector <int> > * screen = new vector <vector <int> >;
    screen->resize(2 * screen_size * screen_resolution);
    for (unsigned int i = 0; i < 2*screen_size*screen_resolution; i++){
        (*screen)[i].resize(2 * screen_size * screen_resolution);
        for (int j = 0; j < 2 * screen_size * screen_resolution; j++){
            (*screen)[i][j] = 0;
        }
    }

    unsigned long long photon_in = 0;
    unsigned long long photon_out = 0;
    int number_of_atoms = 0;
    int number_of_excited_atoms = 0;

    //Массив, заполненный рандомными числами от 0 до 1, нужен для генерирования фотонов от лазера (выбираем рандомную точку в квадрате) 
    // и испускания фотонов от атомов в рандомном направлении (функция emit_photon)
    unsigned long long seed = 5000000; 
    float* random_array = new float[seed];
    unsigned long long random_array_size = seed;
    unsigned long long pointer = 0; 
    mt19937 rng(seed); 
    for (unsigned i = 0; i < random_array_size; i++){
        random_array[i] = float(rng())/(pow(2, 32) - 1);
    }

    std::stringstream stream;
    stream << std::fixed << std::setprecision(2) << intensity_coeff;
    std::string s = stream.str();
    string title = "random_gauss_03.12_56_1_" + s;

    float arr[] = {1., 2., 3.};
    Atom trial_atom(arr, arr, arr, 0, 0, 0, 0, 0, 0.24, false);
    float trial_emission[] = {1, 1, 1};
    float trial_dir[] = {1, 0, 0};
    Photon photon(trial_emission, trial_dir);
    
    // Количество фотонов, которые вылетают из лазера за итерацию
    int photons_per_iter = int(2760 * intensity_coeff * pow(screen_size - 1, 2)/pow(28, 2));

    for (int time = 0; time < 2001; time++){ //TODO поменять время
        std::cout << time << endl;
        //Фотоны вылетают из лазера. 
        for (int i = 0; i < photons_per_iter; i++){ 
            photon_in++;
            if (pointer + 3 > random_array_size){
                pointer = rand()%(random_array_size-5);
            }
            trial_atom.coordinate[0] = -3*screen_size; 
            trial_atom.coordinate[1] = (random_array[pointer++] - 0.5) * 2 * (screen_size-1); 
            trial_atom.coordinate[2] = (random_array[pointer++] - 0.5) * 2 * (screen_size-1); 
            
            photon.dir[0] = 1;
            photon.dir[1] = 0;
            photon.dir[2] = 0;

            std::copy_n(trial_atom.coordinate, 3, photon.emission);
            next_collision(&photon, &trial_atom, 1, ball_cells, ball_cells_size, &photon_out,  screen_size, screen_resolution, screen, screen_len); //TODO убрать trial_atom
        }

        //Проходимся по всем клеткам, в каждой клетке перебираем атомы. Если в прошлой итерации у атома был фотон, излучаем его.
        for (int i = 0; i < ball_cells_size; i++){
            for (int j = 0; j < ball_cells[i].number_of_elements; j++){
                if (ball_cells[i].elements[j].tmp and ball_cells[i].elements[j].photons_inside == 1){
                    ball_cells[i].elements[j].photons_inside = 0;
                    std::copy_n(ball_cells[i].elements[j].coordinate, 3, photon.emission);
                    emit_photon(&photon, &ball_cells[i].elements[j], random_array, random_array_size, &pointer);
                    next_collision(&photon, &ball_cells[i].elements[j], 1, ball_cells, ball_cells_size, &photon_out,  screen_size, screen_resolution, screen, screen_len);
                }
            }
        }

        //Двигаем атомы. Если атом вышел из клетки, помещаем его в новую или в массив атомов, которые вылетели из облака
        for (int i = 0; i < ball_cells_size; i++){
            for (int j = 0; j < ball_cells[i].number_of_elements; j++){
                
                if (ball_cells[i].elements[j].tmp){
                    ball_cells[i].elements[j].move(); 
                    change_cell(&ball_cells[i].elements[j], &ball_cells[i], ball_cells, ball_cells_size, &photon_out, out_of_cloud, size_out_of_cloud, &photon_in);
                }
                if (ball_cells[i].number_of_elements >= cell_size - 2){
                    std::cout<< 100 << endl;
                    exit;
                }
            }
        }
        for (int i = 0; i < *size_out_of_cloud; i++){
            out_of_cloud[i].move();
        }

        for (int i = 0; i < ball_cells_size; i++){
            for (int j = 0; j < ball_cells[i].number_of_elements; j++){
                if (ball_cells[i].elements[j].tmp){
                    ball_cells[i].elements[j].photons_inside = ball_cells[i].elements[j].photons_tmp;
                    ball_cells[i].elements[j].photons_tmp = 0; // TODO segmentation fault
                }
            }
        }

        //Очищаем клетки, которые освободились
        clear_cells(ball_cells, ball_cells_size);

        //Считаем атомы и вылетевшие фотоны для самопроверки
        number_of_atoms = 0;
        number_of_excited_atoms = 0;
        for (int i = 0; i < ball_cells_size; i++){
            for (int j = 0; j < ball_cells[i].number_of_elements; j++){
                if (ball_cells[i].elements[j].tmp == false){
                    continue;
                }
                number_of_atoms++;
                if (ball_cells[i].elements[j].photons_inside == 1) {
                    number_of_excited_atoms++;
                }
                else if (ball_cells[i].elements[j].photons_inside != 0){
                    std::cout << 545 << endl;
                }
            }
        }
        std::cout << intensity_coeff << ' '<< time << ' ' << photon_in << ' ' << photon_out << ' ' << number_of_excited_atoms << ' ' <<  number_of_excited_atoms + photon_out << ' ' << number_of_atoms << ' ' << *size_out_of_cloud << endl;
    
        //Дальше запись в файл  
        if (time % 50 == 0){ 
            write(title + "/movement_data_gauss_"+ to_string(time) + ".txt", ball_cells, ball_cells_size, out_of_cloud, size_out_of_cloud);

            std::ofstream out;
            out.open(title + "/screen_data_gauss_"+ to_string(time) + ".txt");
            if (out.is_open()) {
                for (int i = 0; i < 2 * screen_size * screen_resolution; i++){
                    for (int j = 0; j < 2 * screen_size * screen_resolution; j++){
                        out << (*screen)[i][j] << ' ' << i << ' ' << j << std::endl;
                    }
                }
            }
            out.close();
        } 

    }   
    delete[] random_array;
    delete screen; 

    return;
}

int main(){
    //Массив клеток, в каждой из которого хранится < cell_size атомов
    const int ball_cells_size = 20631; 
    Ball_cell* ball_cells = new Ball_cell[ball_cells_size];

    //Массив вылетевших атомов
    Atom* out_of_cloud = new Atom[10000]; //see change_cell 360
    int size_out_of_cloud = 0;
    
    // Считываем с файла координаты клеток
    ifstream fileptr ("c:\\_current_work\\MIPT\\Programming\\C++\\project\\cells_56.txt"); //Менять размер массива вручную 
    if (not fileptr) cout << "Файл не открыт!!!\n\n"; 
    for (int i = 0; i < ball_cells_size; i++){
        fileptr >> ball_cells[i].coordinate_centre[0] >> ball_cells[i].coordinate_centre[1] >> ball_cells[i].coordinate_centre[2];
        ball_cells[i].edge = 6;
        ball_cells[i].number_of_elements = 0;
    }
    fileptr.close();
    std::cout<< ball_cells[0].coordinate_centre[0]<<endl;

    // Добавляем в клетки атомы
    int number_of_atoms_0 = 73846; //См. cell_size - менять взависимости от концентрации


    Atom trial_atom;
    ifstream file ("c:\\_current_work\\MIPT\\Programming\\C++\\project\\atoms_1_56.txt"); //TODO менять размер массива вручную
   
    if (not file) cout << "файл не открыт!!!\n\n";
    for (int i = 0; i < number_of_atoms_0; i++){ 
        file >> trial_atom.coordinate_0[0] >> trial_atom.coordinate_0[1] >> trial_atom.coordinate_0[2];// >> 
        bool break_flag = true;
        for (int j = 0; j < ball_cells_size; j++){
            if (distance_between_points(trial_atom.coordinate_0, ball_cells[j].coordinate_centre) <  1.72 * ball_cells[j].edge/2){
                std::copy_n(trial_atom.coordinate_0, 3, ball_cells[j].elements[ball_cells[j].number_of_elements].coordinate);
                std::copy_n(trial_atom.coordinate_0, 3, ball_cells[j].elements[ball_cells[j].number_of_elements].coordinate_0);
                ball_cells[j].elements[ball_cells[j].number_of_elements].tmp = true;
                ball_cells[j].number_of_elements++;
                if (ball_cells[j].number_of_elements > cell_size - 2){
                    cout << 469 << endl;
                }
                break_flag = false;
                break;
            }
        }
        if (break_flag){
            std::copy_n(trial_atom.coordinate_0, 3, out_of_cloud[size_out_of_cloud].coordinate);
            std::copy_n(trial_atom.coordinate_0, 3, out_of_cloud[size_out_of_cloud].coordinate_0);
            std::copy_n(trial_atom.impulse, 3, out_of_cloud[size_out_of_cloud].impulse);
            out_of_cloud[size_out_of_cloud].collisions_chaos = trial_atom.collisions_chaos;
            out_of_cloud[size_out_of_cloud].collisions_laser = trial_atom.collisions_laser;
            out_of_cloud[size_out_of_cloud].tmp = true;
            size_out_of_cloud++;
        }

    }
    file.close();

    experiment(7.5, ball_cells, ball_cells_size, out_of_cloud, &size_out_of_cloud);
    return 0;
}