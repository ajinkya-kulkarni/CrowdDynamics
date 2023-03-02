#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <fstream>
#include <numeric>
#include <algorithm>
#include <time.h>

using namespace std;

const float dt = 1e-4;
const int monitor = 500;
float t_final = 410.0;
float t_initial = 0.0;
float time_BELT_starts = 10.0;

// Update in case of self-propelled particles
float LVA = 0.0;
float time_SYS_starts = 0.0;
/**********************************************/

float Cv=20;
const float SPRING=20000;
const float eta = 2;
const int n_ex = 800;
float exp_power;
float exf[n_ex];

const float PI=3.14159;
float gaussStdDevPercRAD = 0.05; // Uncomment this line to specify the standard deviation of the Gaussian w.r.t. mean radius (RAD)

const float SIGMA=213.34;
const float RHO = 100;

const float GRAV=0; //9.81;

float RAD = 0.0025;
float Secondary_RAD = 0.00125;
float SMALL_RAD = 0.001;
float BIG_RAD = 0.004;

const float Wall_RAD = 0.0025;

const float X0 = 0.0;
const float Y0 = 0.0;
const float WIDTH = 0.4;
const float HEIGHT = 0.4;

float r_cutoff = 5.0*(2.0*RAD);
float r_skin = 0.5*RAD;
//float V_surr = PI*pow(r_cutoff,2);

float Xmin = -6.0*Wall_RAD;
float Ymin = -6.0*Wall_RAD;
float Xmax = WIDTH + 6.0*Wall_RAD;
float Ymax = HEIGHT + 6.0*Wall_RAD;

float CONC = 0.00;
const int CORE_PARTICLES = 6084;

int bot_Core_small = 0;
int bot_Core_big = 3042;
int top_Core_small = 0;
int top_Core_big = 3042;
int m = ceil((Ymax - Ymin)/(r_cutoff+r_skin));
int n = ceil((Xmax - Xmin)/(r_cutoff+r_skin));
int ncells = m*n;

vector<vector<int> > nb_part(CORE_PARTICLES);

float BELT_ini_vel = 0.0;

float BELT_vel_x_bot = 0.5;
float BELT_vel_y_right = 0.0;
float BELT_vel_x_top = 0.0;
float BELT_vel_y_left = 0.0;




inline float rand01() {return static_cast<double> (rand())/RAND_MAX;}

inline float randGauss(float sigma)
{
    float phi = rand01()*2*PI;
    float Upsilon = rand01();
    float Psi = -log(Upsilon);
    float r = sigma*sqrt(2.0*Psi);
    float x = r*cos(phi);
    float y = r*sin(phi);
    if (rand01() > 0.5)
        return x;
    else
        return y;
}

/*// for a more accurate calculation of neighbouring particles' volume within the cutoff radius
float circle_int_area(float r1, float r2, float d)
{
    float r = r2;
    float R = r1;
    if(R < r)
    {
        // swap
        r = r1;
        R = r2;
    }
    if (d < R - r)
        return PI*r*r;
    else
    {
        float a1 = r*r*acos((d*d + r*r - R*R)/(2*d*r));
        float a2 = R*R*acos((d*d + R*R - r*r)/(2*d*R));
        float a3 = 0.5*sqrt((-d+r+R)*(d+r-R)*(d-r+R)*(d+r+R));
        return a1 + a2 - a3;
    }
}*/

ofstream ofp_energy("energy.txt", ios::out);
ofstream ofp_monitor("monitor.txt", ios::out);
float wall_energy_input = 0.0;
float sp_pot_energy = 0.0;

class vect
{
public:
    float x, y;
    vect()
        {
        x = y = 0.0f;
        }
    vect(float tx, float ty)
        {
        x = tx;
        y = ty;
        }
    ~vect() {}
    float distance_calc(vect & d2)
    {
        vect d;
        d.x = d2.x - x;
        d.y = d2.y - y;
        return sqrt(d.x * d.x + d.y * d.y);
    }
    void add_vel(const vect & vel_prev, const vect & f, const vect & f_prev, const float m,float dt = 0.0f)
        {
        if (dt == 0.0f)
            {
            x = vel_prev.getX() +  f.getX();
            y = vel_prev.getY() + f.getY();
            }
        else
            {
            x = vel_prev.getX() + ((f.getX()+f_prev.getX())/m)* dt*0.5;
            y = vel_prev.getY() + ((f.getY()+f_prev.getY())/m)* dt*0.5;
            }
        }
    void add_pos_vel(const vect & v, const vect & f, const float m, float dt = 0.0f)
        {
        if (dt == 0.0f)
            {
            x += v.getX();
            y += v.getY();
            }
        else
            {
            x += v.getX() * dt + 0.5*(f.getX()/m)*dt*dt;
            y += v.getY() * dt + 0.5*(f.getY()/m)*dt*dt;
            }
        }
    void add_pos_vel(const vect & v, const float m, float dt = 0.0f)
        {
        if (dt == 0.0f)
            {
            x += v.getX();
            y += v.getY();
            }
        else
            {
            x += (v.getX()/m) * dt;
            y += (v.getY()/m) * dt;
            }
        }

/*    void addc(float f)
        {
        x += f;
        y += f;
        }
    void subtract(const vect & v)
        {
        x -= v.getX();
        y -= v.getY();
        }
    void scale(float mul)
        {
        x *= mul;
        y *= mul;
        }
    void normalize()
        {
        float ln = magnitude();
        x /= ln;
        y /= ln;
        }*/
    float magnitude() const {return sqrt( (x)*(x) + (y)*(y) );}

    float getX() const {return x;}
    float getY() const {return y;}

    void setX(float v) {x = v;}
    void setY(float v) {y = v;}

    };

class particle : public vect
    {
    public:
        vect position, position_prev;
        float displace;
        vect velocity, velocity_prev, vf;
        vect force, force_prev;
        float radius,mass, beta,gamma;
        int species;
        unsigned short wall_tag;
        unsigned short wall_layer;
        float vf_0, v_surr_part, V_surr;
        vect local_order_param;
        float mass_surr = 0, mass_bot = 0, mass_top = 0;
        short unsigned num_surr = 0, num_bot = 0, num_top = 0;


        particle()
            {
            radius = RAD; //creates disks of radius=RAD.
            position.x = rand01()*WIDTH;
            position.y = rand01()*HEIGHT;
            position_prev.x = 0.0;
            position_prev.y = 0.0;
            displace = 0.0;
            mass = SIGMA*PI*radius*radius;
            beta = 0.0;
	    gamma = 0.0;
            velocity.x = 0.0; //initial velocity and force is zero.
            velocity.y = 0.0;
            velocity_prev.x = 0.0; //initial velocity and force is zero.
            velocity_prev.y = 0.0;
            force.x = 0.0;
            force.y = 0.0;
            force_prev.x = 0.0;
            force_prev.y = 0.0;
            v_surr_part = 0;
            vf_0 = 0.0;
            V_surr = PI*pow(5*(2*radius),2);
            species = 0;
            wall_tag = 0;
            wall_layer = 0;
            local_order_param.x = 0.0;
            local_order_param.y = 0.0;
            mass_surr = 0.0; num_surr = 0;
            mass_bot = 0.0; num_bot = 0;
            mass_top = 0.0; num_top = 0;
            }

        particle(float u)
            {
            radius = RAD;
            position.x = rand01()*WIDTH;
            position.y = rand01()*u;
            position_prev.x = 0.0;
            position_prev.y = 0.0;
            displace = 0.0;
            mass = SIGMA*PI*radius*radius;
            beta = 0.0;
	    gamma = 0.0;
            velocity.x = 0.0;
            velocity.y = 0.0;
            velocity_prev.x = 0.0; //initial velocity and force is zero.
            velocity_prev.y = 0.0;
            force.x = 0.0;
            force.y = 0.0;
            force_prev.x = 0.0;
            force_prev.y = 0.0;
            v_surr_part = 0;
            vf_0 = 0.0;
            V_surr = PI*pow(5*(2*radius),2);
            species = 0;            
            wall_tag = 0;
            wall_layer = 0;
            local_order_param.x = 0.0;
            local_order_param.y = 0.0;
            mass_surr = 0.0; num_surr = 0;
            mass_bot = 0.0; num_bot = 0;
            mass_top = 0.0; num_top = 0;
        }

        particle(float u,float v)
            {
            radius = RAD;
            position.x = u;
            position.y = v;
            position_prev.x = u;
            position_prev.y = v;
            displace = 0.0;
            mass = SIGMA*PI*radius*radius;
            beta = 0.0;
	    gamma = 0.0;
            velocity.x = 0.0;
            velocity.y = 0.0;
            velocity_prev.x = 0.0; //initial velocity and force is zero.
            velocity_prev.y = 0.0;
            force.x = 0.0;
            force.y = 0.0;
            force_prev.x = 0.0;
            force_prev.y = 0.0;
            v_surr_part = 0;
            vf_0 = 0.0;
            V_surr = PI*pow(5*(2*radius),2);
            species = 0;
            wall_tag = 0;
            wall_layer = 0;
            local_order_param.x = 0.0;
            local_order_param.y = 0.0;
            mass_surr = 0.0; num_surr = 0;
            mass_bot = 0.0; num_bot = 0;
            mass_top = 0.0; num_top = 0;
        }

        particle(float u,float v, float r)
            {
            radius = r;
            position.x = u;
            position.y = v;
            position_prev.x = u;
            position_prev.y = v;
            displace = 0.0;
            mass = SIGMA*PI*radius*radius;
            beta = 0.0;
	    gamma = 0.0;
            velocity.x = 0.0;
            velocity.y = 0.0;
            velocity_prev.x = 0.0; //initial velocity and force is zero.
            velocity_prev.y = 0.0;
            force.x = 0.0;
            force.y = 0.0;
            force_prev.x = 0.0;
            force_prev.y = 0.0;
            v_surr_part = 0;
            vf_0 = 0.0;
            V_surr = PI*pow(5*(2*radius),2);
            species = 0;
            wall_tag = 0;
            wall_layer = 0;
            local_order_param.x = 0.0;
            local_order_param.y = 0.0;
            mass_surr = 0.0; num_surr = 0;
            mass_bot = 0.0; num_bot = 0;
            mass_top = 0.0; num_top = 0;

        }
        particle(float u,float v, float r, float vel_x, float vel_y)
            {
            radius = r;
            position.x = u;
            position.y = v;
            position_prev.x = u;
            position_prev.y = v;
            displace = 0.0;
            mass = SIGMA*PI*radius*radius;
            beta = 0.0;
	    gamma = 0.0;
            velocity.x = vel_x;
            velocity.y = vel_y;
            velocity_prev.x = 0.0; //initial velocity and force is zero.
            velocity_prev.y = 0.0;
            force.x = 0.0;
            force.y = 0.0;
            force_prev.x = 0.0;
            force_prev.y = 0.0;
            v_surr_part = 0;
            vf_0 = 0.0;
            V_surr = PI*pow(5*(2*radius),2);
            species = 0;
            wall_tag = 0;
            wall_layer = 0;
            local_order_param.x = 0.0;
            local_order_param.y = 0.0;
            mass_surr = 0.0; num_surr = 0;
            mass_bot = 0.0; num_bot = 0;
            mass_top = 0.0; num_top = 0;
        }
	Particle(float u,float v, float r, float vel_x, float vel_y, float p_gamma)
            {
            radius = r;
            position.x = u;
            position.y = v;
            position_prev.x = u;
            position_prev.y = v;
            displace = 0.0;
            mass = SIGMA*PI*radius*radius;
            beta = 0.0;
	    gamma = p_gamma;
            velocity.x = vel_x;
            velocity.y = vel_y;
            velocity_prev.x = 0.0; //initial velocity and force is zero.
            velocity_prev.y = 0.0;
            force.x = 0.0;
            force.y = 0.0;
            force_prev.x = 0.0;
            force_prev.y = 0.0;
            v_surr_part = 0;
            vf_0 = 0.0;
            V_surr = PI*pow(5*(2*radius),2);
            species = 0;
            wall_tag = 0;
            wall_layer = 0;
            local_order_param.x = 0.0;
            local_order_param.y = 0.0;
            mass_surr = 0.0; num_surr = 0;
            mass_bot = 0.0; num_bot = 0;
            mass_top = 0.0; num_top = 0;
        }

        ~particle() {}

        void setPrevPosx(float x_prev) { position_prev.x = x_prev; }
        void setPrevPosy(float y_prev) { position_prev.y = y_prev; }
        void setDisplace(float d) { displace = d; }
        void setVelx(float velx) { velocity.x = velx; }
        void setVely(float vely) { velocity.y = vely; }
        void addCx(float c) { position.x += c; }
        void addCy(float c) { position.y += c; }
        void setRadius(float r) { radius = r; }
        void setMass(float m) { mass = m; }
        void setBeta(float b) { beta = b; }
	void setgamma(float gm) { gamma = gm; }
        void setVsurr(float vsurr) { V_surr = vsurr; }
        void setSpecies(int tag) { species = tag; }
        void setWall(unsigned short tag, unsigned short layer)
        {
            wall_tag = tag;
            wall_layer = layer;
        }

        void collideOther(particle * other, bool VerletUpdate)     //when the disk collides with other disks.
        {
            vect f, dist;
            dist.x = other->position.x - position.x;
            dist.y = other->position.y - position.y;
            float distance = position.distance_calc(other->position);
            float minDist = (other->radius + radius);
            int tmp_id;
            float tmpexf, masstmp;

            if(distance<0.2*minDist)  distance=0.1*minDist;

            // the present particle is in the sphere of influence of the 'other'
            float h = 5*(2*radius);
            if(distance < h)
                {
                tmp_id = 1+n_ex*(distance/h);
                tmpexf = exf[tmp_id];
                masstmp = other->mass*tmpexf;

                v_surr_part += PI*pow(other->radius,2);

                vf.x += masstmp * other->velocity.x;
                vf.y += masstmp * other->velocity.y;
                vf_0 += masstmp;               

                }


            // the 'other' particle is in the sphere of influence of the present
            h = 5*(2*other->radius);
            if (distance < h)
            {
                tmp_id = 1+n_ex*(distance/h);
                tmpexf = exf[tmp_id];
                masstmp = mass*tmpexf;

                other->v_surr_part += PI*pow(radius,2);

                other->vf.x += masstmp * velocity.x;
                other->vf.y += masstmp * velocity.y;
                other->vf_0 += masstmp;

            }


            if (distance < minDist)
            {
                f.x = -(minDist/distance -1) * SPRING * dist.x;
                f.y = -(minDist/distance -1) * SPRING * dist.y;
                if (VerletUpdate)
                    sp_pot_energy += 0.5*SPRING*pow(minDist - distance,2);
            }

            force.x += f.x;
            force.y += f.y;

            other->force.x -= f.x;
            other->force.y -= f.y;
        }
        void orderCalcOther(particle * other)
        {
            vect a, dist;
            float angle;
            dist.x = other->position.x - position.x;
            dist.y = other->position.y - position.y;
            float distance = position.distance_calc(other->position);
            float minDist = (other->radius + radius);

            if (distance < minDist*1.2)
            {
                angle = atan2(dist.y, dist.x);

                local_order_param.x += cos(6*angle);
                local_order_param.y += sin(6*angle);

                angle = atan2(-1.0*dist.y, -1.0*dist.x);
                other->local_order_param.x += cos(6*angle);
                other->local_order_param.y += sin(6*angle);

                num_surr +=1;
                other->num_surr += 1;
                mass_surr += other->mass;
                other->mass_surr += mass;

                if (other->species == 1)
                {
                    num_bot += 1;
                    mass_bot += other->mass;
                }
                else
                {
                    num_top += 1;
                    mass_top += other->mass;
                }

                if (species == 1)
                {
                    other->num_bot += 1;
                    other->mass_bot += mass;
                }
                else
                {
                    other->num_top += 1;
                    other->mass_top += mass;
                }
            }
        }

        void update_pos(const float dt)
            {
            position.add_pos_vel(velocity,force,mass,dt);
            float part_x = position.getX();
            float part_y = position.getY();

            if ( part_x < X0)  {position.setX(RAD*2.0);}
            if ( part_x > X0+WIDTH)  {position.setX(WIDTH - 2.0*RAD);}

            if ( part_y < Y0)  {position.setY(RAD*2.0);}
            if ( part_y > Y0+HEIGHT)  {position.setY(HEIGHT - 2.0*RAD);}
            }

        void update_wall_pos(const float dt)
        {
            position.add_pos_vel(velocity,1.0, dt);
            float part_x = position.getX();
            float part_y = position.getY();

            if (wall_tag == 1 || wall_tag == 3) // Top or bottom wall
            {
                if (wall_layer == 1 || wall_layer == 3) // 1st or 3rd layer
                {
                    if ( part_x > WIDTH + 2*Wall_RAD)  {position.setX(part_x - (WIDTH + 4*Wall_RAD));}
                    if ( part_x < -2*Wall_RAD)  {position.setX(part_x + (WIDTH + 4*Wall_RAD));}
                }
                else                                    // 2nd layer
                {
                    if ( part_x > WIDTH + Wall_RAD)  {position.setX(part_x - (WIDTH + 2*Wall_RAD));}
                    if ( part_x < -Wall_RAD)  {position.setX(part_x + (WIDTH + 2*Wall_RAD));}
                }
            }

            if (wall_tag == 2 || wall_tag == 4) // Side walls (left or right)
            {
                if (wall_layer == 1 || wall_layer == 3) // 1st or 3rd layer
                {
                    if ( part_y > HEIGHT + 2*Wall_RAD)  {position.setY(part_y - (HEIGHT + 4*Wall_RAD));}
                    if ( part_y < -2*Wall_RAD)  {position.setY(part_y + (HEIGHT + 2*Wall_RAD));}
                }
                else                                    // 2nd layer
                {
                    if ( part_y > HEIGHT + Wall_RAD)  {position.setY(part_y - (HEIGHT + 2*Wall_RAD));}
                    if ( part_y < -Wall_RAD)  {position.setY(part_y + (HEIGHT + 2*Wall_RAD));}
                }

            }

        }
        void pseudo_update_vel(const float dt)
        {
            velocity_prev.x = velocity.x;
            velocity_prev.y = velocity.y;
            velocity.add_pos_vel(force, mass, dt);
        }
        void update_vel(const float dt)
        {
            velocity.add_vel(velocity_prev, force, force_prev, mass, dt);
        }

    };

vector<particle> disks;

class simulate : public particle
    {
    public:
        int stop = clock();
        vect gravity;
        vector<int> wall_limits;
        float time = t_initial;
        vect displace_max;
        float power_sum = 0;

        simulate()
            {
            gravity.x = 0.0;
            gravity.y = -GRAV;
            }
        ~simulate() {}

        void boundaryParticlesGenerate()
        {
            ofstream ofp_wall("wall_number.txt", ios::out);
            ofstream ofp_wall_part_stat_1("wall_data_stat_1.txt", ios::out);
            ofstream ofp_wall_part_stat_2("wall_data_stat_2.txt", ios::out);
            ofstream ofp_wall_part_stat_3("wall_data_stat_3.txt", ios::out);
            ofstream ofp_wall_part_stat_4("wall_data_stat_4.txt", ios::out);

            wall_limits.push_back(disks.size());
            //            /***************** Bottom wall (belt) ************************/
            for( float i = -Wall_RAD; i <= WIDTH + Wall_RAD ; i+= RAD * 2 )    // Layer 1
                {
                particle disk = particle(i, -Wall_RAD, Wall_RAD, BELT_ini_vel, 0.0); // Moving wall
                disk.setWall(1,1);    // bottom wall tag = 1
                disks.push_back(disk);
                if (fabs(BELT_vel_x_bot) == 0.0)
                    ofp_wall_part_stat_1 << disk.position.x << "\t" << disk.position.y << "\t"<< disk.radius << "\t"<< disk.velocity.x << "\t" << disk.velocity.y << endl;
                }
            for( float i = 0; i <= WIDTH ; i+= RAD * 2 )    // Layer 2
                {
                particle disk = particle(i, -Wall_RAD*(1+2*sin(PI/3.0)), Wall_RAD, BELT_ini_vel, 0.0); // Moving wall
                disk.setWall(1,2);    // bottom wall tag = 1
                disks.push_back(disk);
                if (fabs(BELT_vel_x_bot) == 0.0)
                    ofp_wall_part_stat_1 << disk.position.x << "\t" << disk.position.y << "\t"<< disk.radius << "\t"<< disk.velocity.x << "\t" << disk.velocity.y << endl;
                }
            for( float i = -Wall_RAD; i <= WIDTH + Wall_RAD ; i+= RAD * 2 )     // Layer 3
                {
                particle disk = particle(i, -Wall_RAD*(1+4*sin(PI/3.0)), Wall_RAD, BELT_ini_vel, 0.0); // Moving wall
                disk.setWall(1,3);    // bottom wall tag = 1
                disks.push_back(disk);
                if (fabs(BELT_vel_x_bot) == 0.0)
                    ofp_wall_part_stat_1 << disk.position.x << "\t" << disk.position.y << "\t"<< disk.radius << "\t"<< disk.velocity.x << "\t" << disk.velocity.y << endl;
                }
            wall_limits.push_back(disks.size());
            unsigned short bot_wall = wall_limits.back() - wall_limits.front();
            ofp_wall << bot_wall << endl;
            ofp_wall_part_stat_1.close();

            //            /***************** Right wall ************************/
            for( float i = HEIGHT + Wall_RAD; i >= -Wall_RAD ; i-= Wall_RAD * 2)     // Layer 1
                {
                particle disk = particle(WIDTH + Wall_RAD, i, Wall_RAD);
                disk.setWall(2,1);    // right wall tag = 2
                disks.push_back(disk);
                if (fabs(BELT_vel_y_right) == 0.0)
                    ofp_wall_part_stat_2 << disk.position.x << "\t" << disk.position.y << "\t"<< disk.radius << "\t"<< disk.velocity.x << "\t" << disk.velocity.y << endl;
                }
            for( float i = HEIGHT; i >= 0.0 ; i-= Wall_RAD * 2)     // Layer 2
                {
                particle disk = particle(WIDTH + Wall_RAD*(1+2*sin(PI/3.0)), i, Wall_RAD);
                disk.setWall(2,2);    // right wall tag = 2
                disks.push_back(disk);
                if (fabs(BELT_vel_y_right) == 0.0)
                    ofp_wall_part_stat_2 << disk.position.x << "\t" << disk.position.y << "\t"<< disk.radius << "\t"<< disk.velocity.x << "\t" << disk.velocity.y << endl;
                }
            for( float i = HEIGHT + Wall_RAD; i >= -Wall_RAD ; i-= Wall_RAD * 2)     // Layer 3
                {
                particle disk = particle(WIDTH + Wall_RAD*(1+4*sin(PI/3.0)), i, Wall_RAD);
                disk.setWall(2,3);    // right wall tag = 2
                disks.push_back(disk);
                if (fabs(BELT_vel_y_right) == 0.0)
                    ofp_wall_part_stat_2 << disk.position.x << "\t" << disk.position.y << "\t"<< disk.radius << "\t"<< disk.velocity.x << "\t" << disk.velocity.y << endl;
                }
            wall_limits.push_back(disks.size());
            unsigned short right_wall = wall_limits.back() - bot_wall;
            ofp_wall << right_wall << endl;
            ofp_wall_part_stat_2.close();


            //            /***************** Top wall ************************/
            for( float i = WIDTH + Wall_RAD; i >= -Wall_RAD ; i-= Wall_RAD * 2)     // Layer 1
                {
                particle disk = particle(i, HEIGHT + Wall_RAD, Wall_RAD);
                disk.setWall(3,1);    // top wall tag = 3
                disks.push_back(disk);
                if (fabs(BELT_vel_x_top) == 0.0)
                    ofp_wall_part_stat_3 << disk.position.x << "\t" << disk.position.y << "\t"<< disk.radius << "\t"<< disk.velocity.x << "\t" << disk.velocity.y << endl;
                }
            for( float i = WIDTH; i >= 0.0 ; i-= Wall_RAD * 2)      // Layer 2
                {
                particle disk = particle(i, HEIGHT + Wall_RAD*(1+2*sin(PI/3.0)), Wall_RAD);
                disk.setWall(3,2);    // top wall tag = 3
                disks.push_back(disk);
                if (fabs(BELT_vel_x_top) == 0.0)
                    ofp_wall_part_stat_3 << disk.position.x << "\t" << disk.position.y << "\t"<< disk.radius << "\t"<< disk.velocity.x << "\t" << disk.velocity.y << endl;
                }
            for( float i = WIDTH + Wall_RAD; i >= -Wall_RAD ; i-= Wall_RAD * 2)     // Layer 3
                {
                particle disk = particle(i, HEIGHT + Wall_RAD*(1+4*sin(PI/3.0)), Wall_RAD);
                disk.setWall(3,3);    // top wall tag = 3
                disks.push_back(disk);
                if (fabs(BELT_vel_x_top) == 0.0)
                    ofp_wall_part_stat_3 << disk.position.x << "\t" << disk.position.y << "\t"<< disk.radius << "\t"<< disk.velocity.x << "\t" << disk.velocity.y << endl;
                }
            wall_limits.push_back(disks.size());
            unsigned short top_wall = wall_limits.back() - bot_wall - right_wall;
            ofp_wall << top_wall << endl;
            ofp_wall_part_stat_3.close();

            //            /***************** Left wall ************************/
            for( float i = -Wall_RAD; i <= HEIGHT + Wall_RAD ; i+= Wall_RAD*2 )     // Layer 1
                {
                particle disk = particle(-Wall_RAD, i, Wall_RAD);
                disk.setWall(4,1);    // left wall tag = 4
                disks.push_back(disk);
                if (fabs(BELT_vel_y_left) == 0.0)
                    ofp_wall_part_stat_4 << disk.position.x << "\t" << disk.position.y << "\t"<< disk.radius << "\t"<< disk.velocity.x << "\t" << disk.velocity.y << endl;
                }
            for( float i = 0.0; i <= HEIGHT ; i+= Wall_RAD*2 )     // Layer 2
                {
                particle disk = particle(-Wall_RAD*(1+2*sin(PI/3.0)), i, Wall_RAD);
                disk.setWall(4,2);    // left wall tag = 4
                disks.push_back(disk);
                if (fabs(BELT_vel_y_left) == 0.0)
                    ofp_wall_part_stat_4 << disk.position.x << "\t" << disk.position.y << "\t"<< disk.radius << "\t"<< disk.velocity.x << "\t" << disk.velocity.y << endl;
                }
            for( float i = -Wall_RAD; i <= HEIGHT + Wall_RAD ; i+= Wall_RAD*2 )     // Layer 3
                {
                particle disk = particle(-Wall_RAD*(1+4*sin(PI/3.0)), i, Wall_RAD);
                disk.setWall(4,3);    // left wall tag = 4
                disks.push_back(disk);
                if (fabs(BELT_vel_y_left) == 0.0)
                    ofp_wall_part_stat_4 << disk.position.x << "\t" << disk.position.y << "\t"<< disk.radius << "\t"<< disk.velocity.x << "\t" << disk.velocity.y << endl;
                }
            wall_limits.push_back(disks.size());
            unsigned short left_wall = wall_limits.back() - bot_wall - right_wall - top_wall;
            ofp_wall << left_wall << endl;
            ofp_wall.close();
            ofp_wall_part_stat_4.close();
        }
        void boundaryParticlesRead()
        {
            int i = 0, wall_tot;
            float x,y,r,velx,vely;
            wall_limits.push_back(disks.size());
            ifstream ifp_wall("wall_ini.txt", ios::in);
            while(!ifp_wall.eof())
            {
                ifp_wall>>x>>y>>r>>velx>>vely;
                particle disk = particle(x,y,r,velx,vely);
                disks.push_back(disk);
                i++;
            }                        
            disks.pop_back();
            wall_tot = disks.size();
            // Assigning wall tags in a single statement
            //for (int j=0; j <disks.size(); j++)
              //  disks[j].setWall((j+1)%(disks.size()/4),1);    // The read particles are assigned tags 1, 2, 3, 4

            // The read particles are assigned tags 1, 2, 3, 4
            // Assigning wall tags in 4 statements
            for (int j=0; j <disks.size()/4; j++)
                disks[j].setWall(1,1);
            for (int j=disks.size()/4; j <disks.size()/2; j++)
                disks[j].setWall(2,1);
            for (int j=disks.size()/2; j <disks.size()*3/4; j++)
                disks[j].setWall(3,1);
            for (int j=disks.size()*3/4; j <disks.size(); j++)
                disks[j].setWall(4,1);

            wall_limits.push_back(disks.size()/4);
            wall_limits.push_back(disks.size()/2);
            wall_limits.push_back(disks.size()*3/4);
            wall_limits.push_back(disks.size());
        }

        void interiorParticlesGenerate(int i)
        {            
            // N_s = N* (x / (1-x*(r/R)^2))
            int small_part = floor(CORE_PARTICLES*CONC/(1-CONC*(1-pow(Secondary_RAD/RAD,2))));

            // N_b = N* ( (1-x) / ( 1 - x * (1 - (r/R)^2) ) )
            int big_part = floor(CORE_PARTICLES*(1 - CONC)/(1-CONC*(1-pow(Secondary_RAD/RAD,2))));
            ofstream ofp_int("bot_top_small_big_number.txt", ios::out);

            bot_Core_small = floor(small_part/2);
            bot_Core_big = floor(big_part/2);
            top_Core_small = small_part - bot_Core_small;
            top_Core_big = big_part - bot_Core_big;

            for (int i = 0; i < bot_Core_small; i++)
                {
                particle disk = particle(HEIGHT/2.0);
                disk.setSpecies(1);
                disk.setRadius(Secondary_RAD);
                disk.setMass(SIGMA*PI*pow(Secondary_RAD,2));
                disk.setVsurr(PI*pow(5*2*Secondary_RAD,2));
                disks.push_back(disk);
                }
            ofp_int << bot_Core_small << endl;
            for (int i = bot_Core_small; i < bot_Core_small+bot_Core_big; i++)
                {
                particle disk = particle(HEIGHT/2.0);
                disk.setSpecies(1);
                disks.push_back(disk);
                }
            ofp_int << bot_Core_big << endl;
            for (int i = bot_Core_small+bot_Core_big; i < bot_Core_small+bot_Core_big + top_Core_small; i++)
                {
                particle disk = particle(HEIGHT/2.0);
                disk.addCy(HEIGHT/2.0);
                disk.setSpecies(2);
                disk.setRadius(Secondary_RAD);
                disk.setMass(SIGMA*PI*pow(Secondary_RAD,2));
                disk.setVsurr(PI*pow(5*2*Secondary_RAD,2));
                disks.push_back(disk);
                }
            ofp_int << top_Core_small << endl;
            for (int i = bot_Core_small+bot_Core_big + top_Core_small; i < bot_Core_small+bot_Core_big + top_Core_small + top_Core_big; i++)
                {
                particle disk = particle(HEIGHT/2.0);
                disk.addCy(HEIGHT/2.0);
                disk.setSpecies(2);
                disks.push_back(disk);
                }
            ofp_int << top_Core_big << endl;
            ofp_int.close();
            if (i>1)
            {
                float norm_area = CORE_PARTICLES * PI*RAD*RAD, tot_area;
                vector<float> r;
                if (i == 2)           // Uniform random
                {
                for (int i = 0; i< CORE_PARTICLES; i++)
                    r.push_back(rand01()*(BIG_RAD - SMALL_RAD) + SMALL_RAD);
                }

                if (i == 3)// Gaussian random
                {

                    float norm_area = CORE_PARTICLES * PI*RAD*RAD, tot_area;
                            // Gaussian random
                    float min = 0.0, max = 0.0;
                    for (int i = 0; i< CORE_PARTICLES; i++)
                    {
                        r.push_back(randGauss(1.0));
                        if (r.back() < min) min = r.back();
                        if (r.back() > max) max = r.back();
                    }
                    float mean = accumulate(r.begin(), r.end(), 0.0)/r.size();
                    float std_dev = 0.0;
                    for_each (r.begin(), r.end(), [&](const float rad) { std_dev += (rad - mean) * (rad - mean); });
                    std_dev = sqrt(std_dev/(r.size() - 1));

                    for (int i = 0; i< CORE_PARTICLES; i++)
                        r[i] = (r[i] - mean) / (max - min); // Mean normalization and scaling done

                    mean = accumulate(r.begin(), r.end(), 0.0)/r.size();
                    std_dev = 0.0;
                    for_each (r.begin(), r.end(), [&](const float rad) { std_dev += (rad - mean) * (rad - mean); });
                    std_dev = sqrt(std_dev/(r.size() - 1));

                    for (int i = 0; i< CORE_PARTICLES; i++)
                        r[i] = r[i]*gaussStdDevPercRAD * RAD/std_dev + RAD;
                }
                tot_area = accumulate(r.begin(), r.end(), 0.0, [](float total, float rad){return total + PI*rad*rad;});

                float max_RAD = RAD;
                for (int i = 0; i< CORE_PARTICLES; i++)
                {
                    r[i] *= sqrt(norm_area/tot_area);
                    disks[i+wall_limits.back()].setRadius(r[i]);
                    disks[i+wall_limits.back()].setMass(SIGMA*PI*r[i]*r[i]);
                    disks[i+wall_limits.back()].setVsurr(PI*pow(5*2*r[i],2));
                    if (r[i] > max_RAD)
                        max_RAD = r[i];
                }
                r_cutoff = 5.0*(2.0*max_RAD);
                r_skin = 0.5*RAD;
            }
        }

       /* void interiorParticlesRead()
        {
            int i = 0;
            float x,y,r,velx,vely, lva;
            ifstream ifp_int("interior_ini.txt", ios::in);
            while(!ifp_int.eof())
            {
                ifp_int>>x>>y>>r>>velx>>vely>>lva;
                particle disk = particle(x,y,r,velx,vely);
                if (i < bot_Core_small + bot_Core_big)
                    disk.setSpecies(1);
                else
                    disk.setSpecies(2);
                disks.push_back(disk);
                i++;
            }
            disks.pop_back();
        }*/
	
	void interiorParticlesRead()
        {
            int i = 0;
            float x,y,r,velx,vely, p_gamma;
            ifstream ifp_int("interior_ini.txt", ios::in);
            while(!ifp_int.eof())
            {
                ifp_int>>x>>y>>r>>velx>>vely>>p_gamma;
                Particle disk = Particle(x,y,r,velx,vely, p_gamma);
                if (i < bot_Core_small + bot_Core_big)
                    disk.setSpecies(1);
                else
                    disk.setSpecies(2);
                disks.push_back(disk);
                i++;
            }
            disks.pop_back();
        }

        void neighbour_list()
        {            
            nb_part.resize(disks.size());
            vector<vector<int> > cell_part(ncells);

        // Un-comment the below 4 lines in case the domain is expanding
            n = ceil((Xmax - Xmin)/(r_cutoff+r_skin)); // n has been updated
            m = ceil((Ymax - Ymin)/(r_cutoff+r_skin)); // m has been updated
            ncells = n*m; // ncells has been updated
            vector<int> head(ncells);
            vector<int> tail(disks.size());

            int cell_ix, cell_iy, cell_n;
            int l;

            for (unsigned short j=0; j < ncells; j++)  head[j] = -1;

            // Head-tail construction
            for (int i = 0 ; i < disks.size(); i++)
                {
                nb_part[i].clear();
                cell_ix = floor((disks[i].position.x - Xmin)/(r_cutoff+r_skin));
                cell_iy = floor((disks[i].position.y - Ymin)/(r_cutoff+r_skin));
                cell_n = cell_ix + cell_iy * n;
                tail[i] = head[cell_n];
                head[cell_n] = i;

                disks[i].setPrevPosx(disks[i].position.x);
                disks[i].setPrevPosy(disks[i].position.y);
                }

            // particles in cell construction
            for (int i=0; i < ncells; i++)
            {
                if (head[i] > -1)
                {
                    cell_part[i].clear(); //cell_part.shrink_to_fit();

                    l=0;
                    cell_part[i].push_back(head[i]);
                    while (tail[cell_part[i][l]] > -1)
                        {
                        cell_part[i].push_back(tail[cell_part[i][l]]);
                        l++;
                        }
                    /*
                      cell_part[i] contains the ids of particles in the i-th cell.
                      The interior particles have higher ids and are more probable to be heads in a given cell.
                    */
                }
            }


            // Neigbourhood list construction
            int potential_nb_cell;
            vector <int> nb_cells;
            vector<int> nb_add(4);
            nb_add[0] = n; nb_add[1] = n+1; nb_add[2] = 1; nb_add[3] = -n+1;

            for (int i=0; i < ncells; i++)
            {
                if (head[i] > -1)
                {

                    for (int j = 0; j< cell_part[i].size(); j++)
                    {
                        for (int k = j+1; k < cell_part[i].size(); k++)
                        {
                            if (disks[cell_part[i][j]].wall_tag == 0) // In order to avoid checking for interaction between wall particles.
                            {
                                if (disks[cell_part[i][j]].position.distance_calc(disks[cell_part[i][k]].position) < r_cutoff + r_skin)
                                    nb_part[cell_part[i][j]].push_back(cell_part[i][k]);
                            }
                        }
                    }
                    /************************************************************************

                        The id of wall_particles in the domain is lower than that of interior ones.
                        Therefore, if cell_part[i][j] is the id of a wall particle,
                        it should not check for interaction with the subsequent wall particles
                    ************************************************************************/

                    nb_cells.clear();
                    for (auto potential_nb_cell_add: nb_add)
                    {
                        potential_nb_cell = potential_nb_cell_add + i;
                        if ((potential_nb_cell < ncells) && (potential_nb_cell > -1))
                        {
                            if (i%n == 0)
                            {
                                nb_cells.push_back(potential_nb_cell);
                            }
                            else
                            {
                                if (potential_nb_cell % n != 0)
                                    nb_cells.push_back(potential_nb_cell);
                            }
                        }
                    }

                    for (auto part: cell_part[i])
                    {
                        for (auto nb_cell: nb_cells)
                        {
                            for (auto neighbour_particle: cell_part[nb_cell])
                            {
                                if (!(disks[part].wall_tag != 0 && disks[neighbour_particle].wall_tag != 0))
                                    // In order to avoid checking for interaction between wall particles.
                                    // We require nb_part[part] to be populated basically; but don't this (populate nb_part[part]) only if both
                                    // the particle in the cell (part) and the neighbouring particle
                                    // belong to the wall.
                                {
                                    float r_cut_loc = 5*2*disks[part].radius;
                                    if (r_cut_loc < SMALL_RAD + BIG_RAD) r_cut_loc = SMALL_RAD + BIG_RAD; // updating local r_cutoff in order to at least calculate the spring force
                                    if (disks[part].position.distance_calc(disks[neighbour_particle].position) < r_cut_loc + r_skin)
                                        // The above conditional is for the current implementation of fluid drag force (5 times particle diameter)
                                        // In case, a common r_cutoff is to be used for drag calculation for all particles,
                                        // comment the above conditional and uncomment the below one.
                                    //if (disks[part].position.distance_calc(disks[neighbour_particle].position) < r_cutoff + r_skin)
                                        nb_part[part].push_back(neighbour_particle);
                                }
                            }
                        }

                    }
                }
            }
        }
        void order_calc_init()
        {
            for(unsigned short i = wall_limits.back() ; i < disks.size(); i++)
            {
                disks[i].num_top = 0;
                disks[i].num_bot = 0;
                disks[i].num_surr = 0;
                disks[i].mass_top = 0.0;
                disks[i].mass_bot = 0.0;
                disks[i].mass_surr = 0.0;
                disks[i].local_order_param.x = 0.0;
                disks[i].local_order_param.y = 0.0;
            }
        }
        void order_calc()
        {
            for(int i = 0; i<disks.size(); i++)
                for (auto nb: nb_part[i])
                    disks[i].orderCalcOther(&disks[nb]);
        }

        void force_initialize(bool VerletUpdate)
        {
            if (VerletUpdate) sp_pot_energy = 0.0;
            if (!VerletUpdate) displace_max.x = 0.0; // initialize the maximum displacement to 0 at the second step of Velocity Verlet

            for(unsigned short i = 0 ; i < disks.size(); i++)
            {
                disks[i].force.x = disks[i].mass*gravity.x;
                disks[i].force.y = disks[i].mass*gravity.y;
                disks[i].vf.x = 0.0;
                disks[i].vf.y = 0.0;
                disks[i].vf_0 = 0.0;
                disks[i].v_surr_part = 0.0;
            }
        }

        void force_calc(bool VerletUpdate)
        {
            for(int i = 0; i<disks.size(); i++)
            {
                for (auto nb: nb_part[i])
                    disks[i].collideOther(&disks[nb], VerletUpdate);

                /********Additional step to calculate the maximum displacement among all the particles*********/
                if (!VerletUpdate) // find the maximum displacement only at the second step of Velocity Verlet
                    if (disks[i].displace > displace_max.x)
                    {
                        displace_max.x = disks[i].displace;
                        displace_max.y = i;
                    }


                /**********************************************************************************************/
            }
        }

        void force_update(bool VerletUpdate)
        {
            float mfrac = 1.0, m_fluid, m_surr, vmag;
            for(unsigned short i = wall_limits.back(); i < disks.size(); i++)
            {
                m_fluid = RHO * (disks[i].V_surr - disks[i].v_surr_part);

                m_surr = SIGMA * disks[i].v_surr_part;
                if (m_fluid<0)
                {
                    m_fluid=0.0;
                }
                mfrac = m_surr/(m_surr + m_fluid);

                if (disks[i].vf_0 > 0)
                {
                    disks[i].vf.x /= disks[i].vf_0;
                    disks[i].vf.y /= disks[i].vf_0;
                }
                else
                {
                    disks[i].vf.x = 0;
                    disks[i].vf.y = 0;
                }

                vmag = disks[i].velocity.magnitude();
                if (vmag == 0.0) vmag =1.0;
                disks[i].force.x += Cv*(2*disks[i].radius)*(mfrac * disks[i].vf.x - disks[i].velocity.x) + disks[i].mass*disks[i].beta*disks[i].velocity.x/vmag - disks[i].mass*disks[i].gamma*(sin(theta_gamma));
                disks[i].force.y += Cv*(2*disks[i].radius)*(mfrac * disks[i].vf.y - disks[i].velocity.y) + disks[i].mass*disks[i].beta*disks[i].velocity.y/vmag + disks[i].mass*disks[i].gamma*(cos(theta_gamma));

                if (VerletUpdate)
                {
                    disks[i].force_prev.x = disks[i].force.x;
                    disks[i].force_prev.y = disks[i].force.y;
                }

            }
        }

        void render() // numerical integration (Velocity-Verlet algorithm)
        {
            force_initialize(true);
            force_calc(true);
            force_update(true);

            for(unsigned short i = wall_limits.back(); i < disks.size(); i++)
            {
                disks[i].update_pos(dt);
                disks[i].pseudo_update_vel(dt);

                disks[i].setDisplace(disks[i].position.distance_calc(disks[i].position_prev));
                //cout << "Displacement = " << disks[i].displace<< endl;
            }

            power_sum = 0;
            for(unsigned short i = wall_limits.front(); i < wall_limits.back(); i++)
            {
                if (fabs(disks[i].velocity.x) > 0 || fabs(disks[i].velocity.y) > 0)
                    disks[i].update_wall_pos(dt);
                power_sum += (disks[i].velocity.x * disks[i].force.x + disks[i].velocity.y * disks[i].force.y);
                disks[i].setDisplace(disks[i].position.distance_calc(disks[i].position_prev));

            }


            force_initialize(false);
            force_calc(false);
            force_update(false);

            for(unsigned short i = wall_limits.back(); i < disks.size(); i++)
            {
                disks[i].update_vel(dt);
            }

        }

        void mainLoop()
        {
        char time_data[50], time_data_1[50],time_data_2[50],time_data_3[50],time_data_4[50];
        int start =0;
        short unsigned nb_list_calc = 0;

        double duration, mix_num = 0.0, mix_mass = 0.0, global_order_param = 0.0, local_param_mag = 0.0;
        int num_iterations = ceil((t_final - t_initial)/dt);
        cout << "Total iterations = " << num_iterations << endl;

        neighbour_list();
        force_initialize(true);
        force_calc(true);
        order_calc_init();
        order_calc();
        nb_list_calc++;

        for (int j = 0; j <= num_iterations; j++)
        {
        	
        	
        	
            wall_energy_input += (-power_sum)*dt;

            if (displace_max.x > 0.99*r_skin) // Update the neighbour list only if the maximum displacement is close to the "skin" thickness for the VERLET list
            {                
                //cout << "Max diplacement = " << displace_max.x <<" at " << displace_max.y << " Wall tag = " << disks[displace_max.y].wall_tag<< endl;
                neighbour_list();
                nb_list_calc++;
            }

            if ((LVA>0) && ((time > time_SYS_starts && time < time_SYS_starts + 2*dt) || (t_initial > time_SYS_starts && time < t_initial + 2*dt)))
            {
                float beta_particle = 1.0;
                Cv = beta_particle*LVA; // used LVA as cv because beta = 1

                for(unsigned short i = wall_limits.back(); i < disks.size(); i++)
                {
                    disks[i].setBeta(beta_particle);
                }
                cout << "Self_propulsion starts with BETA = " << beta_particle << endl;
                cout << "Coordination set to = " << Cv << endl << endl;
            }

            if ((time > time_BELT_starts && time < time_BELT_starts + 2*dt) || (t_initial > time_BELT_starts && time < t_initial + 2*dt))
            {

                for(unsigned short i = wall_limits.front(); i < wall_limits.back(); i++)
                {
                	 
                    if (disks[i].wall_tag == 1) // this is redundant if only one wall is moving
                    {
                    	                  
                        disks[i].setVelx(BELT_vel_x_bot);
                    }

                    if (disks[i].wall_tag == 2) // this is redundant if only one wall is moving
                        disks[i].setVely(BELT_vel_y_right);

                    if (disks[i].wall_tag == 3) // this is redundant if only one wall is moving
                        disks[i].setVelx(BELT_vel_x_top);

                    if (disks[i].wall_tag == 4) // this is redundant if only one wall is moving
                        disks[i].setVely(BELT_vel_y_left);

                }
                cout << "Belt starts " << endl;
            }

            if (j % monitor == 0 )
            {
                float grav_potential_energy =0.0, tot_kinetic_energy = 0.0;
                
                
                snprintf(time_data,sizeof time_data, "./particle_data/interior_%.2f.txt",time);
                ofstream ofp1(time_data, ios::out); // ofp is the identifier. ios::in means it is for writing
                order_calc_init();
                order_calc();
                for(unsigned short i = wall_limits.back(); i < disks.size(); i++)
                {
                    grav_potential_energy += disks[i].mass * GRAV * disks[i].position.y;
                    tot_kinetic_energy += 0.5 * disks[i].mass * pow(disks[i].velocity.magnitude(),2);

                    if (disks[i].num_surr > 0)
                    {
                        disks[i].local_order_param.x /= disks[i].num_surr;
                        disks[i].local_order_param.y /= disks[i].num_surr;
                        mix_num += pow( (disks[i].num_top - disks[i].num_bot)/disks[i].num_surr , 2);
                    }
                    if (disks[i].mass_surr > 0.0)
                    {
                        mix_mass += pow( (disks[i].mass_top - disks[i].mass_bot)/disks[i].mass_surr , 2);
                    }
                    local_param_mag = disks[i].local_order_param.magnitude();
                    ofp1 << disks[i].position.x << "\t" << disks[i].position.y << "\t"<< disks[i].radius << "\t"<< disks[i].velocity.x << "\t" << disks[i].velocity.y << "\t"<< local_param_mag <<endl;
                    global_order_param += local_param_mag;
                }
                ofp1.close();
                mix_num = sqrt(mix_num/(disks.size() - wall_limits.back()));
                mix_mass = sqrt(mix_mass/(disks.size() - wall_limits.back()));
                global_order_param /= disks.size() - wall_limits.back();

                ofp_monitor<<time<<"\t"<<mix_num<<"\t"<<mix_mass<<"\t"<<global_order_param<<endl;

                snprintf(time_data_1,sizeof time_data_1, "./particle_data/wall_1_%.2f.txt",time);
                snprintf(time_data_2,sizeof time_data_2, "./particle_data/wall_2_%.2f.txt",time);
                snprintf(time_data_3,sizeof time_data_3, "./particle_data/wall_3_%.2f.txt",time);
                snprintf(time_data_4,sizeof time_data_4, "./particle_data/wall_4_%.2f.txt",time);

                ofstream ofp2_1(time_data_1, ios::out);
                ofstream ofp2_2(time_data_2, ios::out);
                ofstream ofp2_3(time_data_3, ios::out);
                ofstream ofp2_4(time_data_4, ios::out);



                for(unsigned short i = wall_limits.front(); i < wall_limits.back(); i++)
                {
			
			
                    if (BELT_vel_x_bot != 0)
                    {
                        if (disks[i].wall_tag == 1)
                            ofp2_1 << disks[i].position.x << "\t" << disks[i].position.y << "\t"<< disks[i].radius << "\t"<< disks[i].velocity.x << "\t" << disks[i].velocity.y << endl;
                    }
                    if (BELT_vel_y_right != 0)
                    {
                        if (disks[i].wall_tag == 2)
                            ofp2_2 << disks[i].position.x << "\t" << disks[i].position.y << "\t"<< disks[i].radius << "\t"<< disks[i].velocity.x << "\t" << disks[i].velocity.y << endl;
                    }
                    if (BELT_vel_x_top != 0)
                    {
                        if (disks[i].wall_tag == 3)
                            ofp2_3 << disks[i].position.x << "\t" << disks[i].position.y << "\t"<< disks[i].radius << "\t"<< disks[i].velocity.x << "\t" << disks[i].velocity.y << endl;
                    }
                    if (BELT_vel_y_left != 0)
                    {
                        if (disks[i].wall_tag == 4)
                            ofp2_4 << disks[i].position.x << "\t" << disks[i].position.y << "\t"<< disks[i].radius << "\t"<< disks[i].velocity.x << "\t" << disks[i].velocity.y << endl;
                    }
                }
                ofp2_1.close();
                ofp2_2.close();
                ofp2_3.close();
                ofp2_4.close();

                ofp_energy << time << "\t" << -power_sum << "\t" << wall_energy_input << "\t" << grav_potential_energy << "\t" << tot_kinetic_energy << "\t" << sp_pot_energy << "\n";


                start=stop;
                stop=clock();
                duration = (double)(stop - start) / CLOCKS_PER_SEC;
                cout << "Time = " << time << "\t" << "Simulation time = " << duration << endl;
                cout << "Total no. of particles = " << disks.size() << "\t wall = " << wall_limits.back() << "\t interior = " << disks.size() - wall_limits.back() << endl;
                cout << "Total Verlet list updates = " << nb_list_calc << endl<< endl;
                nb_list_calc = 0;

            }            

            time = (j+1)*dt + t_initial;
            // modify for oscilation
		
		 if ((time > time_BELT_starts )&& (j % monitor == 0 ))
            	{
            		//BELT_vel_x_bot= 0.5 - 1*sin(time/2); // modify the belt velocity
		        for(unsigned short i = wall_limits.front(); i < wall_limits.back(); i++)
		        {
		        	 
		            if (disks[i].wall_tag == 1) // this is redundant if only one wall is moving
		            {
		            	                  
		                disks[i].setVelx(BELT_vel_x_bot);
		            }

		            if (disks[i].wall_tag == 2) // this is redundant if only one wall is moving
		                disks[i].setVely(BELT_vel_y_right);

		            if (disks[i].wall_tag == 3) // this is redundant if only one wall is moving
		                disks[i].setVelx(BELT_vel_x_top);

		            if (disks[i].wall_tag == 4) // this is redundant if only one wall is moving
		                disks[i].setVely(BELT_vel_y_left);

		        }
		        cout << "Belt velocity = " << BELT_vel_x_bot << endl;
            	}
		
            render();
         }
        }
    };



int main (int argc,char *argv[])
{
    bool interiorParticlesProvided = 0, wallParticlesProvided = 0;
    short unsigned interiorPartType = 1;
    if( argc > 1 )
    {
      t_initial = atof(argv[1]);
      cout << "\n ***** Initial time = " << t_initial << " s"<< endl;
      t_final = atof(argv[2]);
      cout << "\n ***** Simulation until time = " << t_final << " s" << endl;
      interiorParticlesProvided = atof(argv[3]);
      wallParticlesProvided = atof(argv[4]);

      RAD = atof(argv[5]); // Not valid if initial interior particles' file is provided
      CONC = atof(argv[6]); // Not valid if initial interior particles' file is provided
      Secondary_RAD = atof(argv[7]); // Not valid if initial interior particles' file is provided

      interiorPartType = atof(argv[8]); // Not valid if initial interior particles' file is provided. Default value of 1 corresponds to mono-disperse or concentration based distribution
      gaussStdDevPercRAD = atof(argv[9]); // Not valid if initial interior particles' file is provided

      SMALL_RAD = atof(argv[10]); // Not valid if initial interior and wall particles' file is provided
      BIG_RAD = atof(argv[11]); // Not valid if initial interior and wall particles' file is provided

      BELT_vel_x_bot = atof(argv[12]);
      BELT_vel_y_left = atof(argv[13]);
      BELT_vel_x_top = atof(argv[14]);
      BELT_vel_y_right = atof(argv[15]);
      time_BELT_starts = atof(argv[16]);
      cout << "\n ***** Belt starts at = " << time_BELT_starts << " s" << endl;

      LVA = atof(argv[17]);
      time_SYS_starts = atof(argv[18]);
      cout << "\n ***** Local Velocity Anarchy = " << LVA << endl;
    }


    srand(unsigned (time(0))); //initialize randomization

    // Generating look-up table
    for (unsigned int j=0; j <= n_ex; j++)
    {
        exp_power = static_cast<float> (j)/n_ex;
        exf[j] = exp(-eta * exp_power * exp_power);
    }

    simulate start = simulate();

    if (wallParticlesProvided)        start.boundaryParticlesRead();
    else        start.boundaryParticlesGenerate();


    if (interiorParticlesProvided)    start.interiorParticlesRead();
    else        start.interiorParticlesGenerate(interiorPartType);

    start.mainLoop();

    ofp_energy.close();
    ofp_monitor.close();
    return 0;
}
