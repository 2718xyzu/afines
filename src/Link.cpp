/*
 *  Link.cpp
 *  
 *
 *  Created by Simon Freedman on 9/10/2014
 *  Copyright 2014 University of Chicago. All rights reserved.
 *
 */

#include "Link.h"
#include "globals.h"
#include "filament.h"

Link::Link(){ }

Link::Link(double len, double stretching_stiffness, filament* f, 
        array<int, 2> myaindex, array<double, 2> myfov, array<int, 2> mynq)
{
    kl      = stretching_stiffness;
    l0      = len;
    fil     = f;
    aindex  = myaindex;
    fov     = myfov;
    nq      = mynq;

    hx = {0,0};
    hy = {0,0};

    force = 0;
    this->step();
}
Link::~Link(){ 
    //std::cout<<"DELETING LINK\n";
};

array<double,2> Link::get_hx(){
    return hx;
}

array<double,2> Link::get_hy(){
    return hy;
}

void Link::set_aindex1(int i){
    aindex[1] = i;
    this->step();
}

// stepping kinetics
void Link::step()
{
    hx[0] = fil->get_actin(aindex[0])->get_xcm();
    hx[1] = fil->get_actin(aindex[1])->get_xcm();
    hy[0] = fil->get_actin(aindex[0])->get_ycm();
    hy[1] = fil->get_actin(aindex[1])->get_ycm();

    xcm = (hx[0]+hx[1])/2.0;
    ycm = (hy[0]+hy[1])/2.0;
    phi=atan2(hy[1]-hy[0],hx[1]-hx[0]);

}

void Link::update_force(string bc, double shear_dist)
{
    force = kl * (dist_bc(bc, hx[1]-hx[0], hy[1]-hy[0], fov[0], fov[1], shear_dist) - l0);
}

double Link::get_force()
{
    return force;
}

void Link::filament_update()
{
    double fx0, fy0, fx1, fy1;
    
    fx0 =  force * cos(phi); 
    fy0 =  force * sin(phi); 
    fx1 =  -fx0;
    fy1 =  -fy0;
    fil->update_forces(aindex[0], fx0, fy0);
    fil->update_forces(aindex[1], fx1, fy1);

}

double Link::get_kl(){
    return kl;
}

double Link::get_l0(){
    return l0;
}
double Link::get_xcm(){
    return xcm;
}

double Link::get_ycm(){
    return ycm;
}

double Link::get_angle(){
    return phi;
}

double Link::get_length(){
    return l0; 
}

std::string Link::write(){
    return "\n" + std::to_string(hx[0]) + "\t" + std::to_string(hy[0]) + "\t" + std::to_string(hx[1]-hx[0]) + "\t" 
        + std::to_string(hy[1]-hy[0]);
}

std::string Link::to_string(){
    
    char buffer [100];
    sprintf(buffer, "aindex[0] = %d;\t aindex[1] = %d;\t kl = %f;\t l0 = %f\nfilament : \n",
                        aindex[0], aindex[1], kl, l0);

    return buffer + fil->to_string();

}

bool Link::operator==(const Link& that) 
{
    /*Note: you can't compare the filament objects because that will lead to infinite recursion;
     * this function requires the filament poiner to be identical to evaluate to true*/
    return (this->aindex[0] == that.aindex[0] && this->aindex[1] == that.aindex[1] &&
            this->kl == that.kl && 
            this->l0 == that.l0 && this->fil == that.fil);
}

bool Link::is_similar(const Link& that) 
{
    
    /* Same as ==; but doesn't compare the filament pointer*/

    return (this->aindex[0] == that.aindex[0] && this->aindex[1] == that.aindex[1] &&
            this->kl == that.kl &&
            this->l0 == that.l0);
}

//All these things swtich from being applicable to actin to being applicable to links:
// Updates all derived quantities of a monomer
void Link::quad_update(string bc, double delrx){
    
    //quadrant numbers crossed by the actin in x-direction
    quad.clear();
    int xlower, xupper, ylower, yupper;
    
    if(hx[0] <= hx[1])
    {
        xlower = int(round(hx[0]/fov[0]*nq[0]));
        xupper = max(int(round( hx[1]/fov[0]*nq[0])), xlower + 1);
    }
    else
    {
        xlower = int(round(hx[1]/fov[0]*nq[0]));
        xupper = max(int(round( hx[0]/fov[0]*nq[0])), xlower + 1);
    };
    
    if(hy[0] <= hy[1])
    {
        ylower = int(round(hy[0]/fov[1]*nq[1]));
        yupper = max(int(round( hy[1]/fov[1]*nq[1])), ylower + 1);
    }
    else
    {
        ylower = int(round(hy[1]/fov[1]*nq[1]));
        yupper = max(int(round( hy[0]/fov[1]*nq[1])), ylower + 1);
    };
    for(int xcoord = xlower; xcoord <= xupper; xcoord++)
    {
        for(int ycoord = ylower; ycoord <= yupper; ycoord++)
        {
            quad.push_back({xcoord, ycoord});
            
            if (abs(xcoord) == nq[0]/2 || abs(ycoord) == nq[1]/2){
                if (bc == "PERIODIC"){
                    if (abs(ycoord) !=  nq[1]/2)                // at xboundary 
                        quad.push_back({-xcoord,  ycoord});
                    else if (abs(xcoord) !=  nq[0]/2)           // at yboundary
                        quad.push_back({ xcoord, -ycoord});
                    else{                                       // at corner
                        quad.push_back({-xcoord,  ycoord});
                        quad.push_back({ xcoord, -ycoord});
                        quad.push_back({-xcoord, -ycoord});
                    }
                }
                else if (bc == "XPERIODIC"){
                    if (abs(xcoord) ==  nq[0]/2) 
                        quad.push_back({-xcoord,  ycoord});
                }
                else if (bc == "LEES-EDWARDS"){
                    if (abs(ycoord) !=  nq[1]/2) 
                        quad.push_back({-xcoord,  ycoord});
                    else if (abs(xcoord) !=  nq[0]/2) 
                        quad.push_back({ xcoord - sgn(ycoord)*int(round(delrx)), -ycoord});
                    else{
                        quad.push_back({-xcoord,  ycoord});
                        quad.push_back({ xcoord - sgn(ycoord)*int(round(delrx)), -ycoord});
                        if (ycoord > 0){
                            if (xcoord - int(round(delrx)) < nq[0]/2)
                                quad.push_back( {xcoord - int(round(delrx)) + 1, -ycoord} );
                            else
                                quad.push_back( {-xcoord, -ycoord} );
                        }
                        else{
                            if (xcoord + int(round(delrx)) > -nq[0]/2)
                                quad.push_back( {xcoord + int(round(delrx)) - 1, -ycoord} );
                            else
                                quad.push_back( {-xcoord, -ycoord} );
                        }
                    }
                }
                     
            }
        }
    }
}

//shortest(perpendicular) distance between an arbitray point and the Link
double Link::get_distance(string bc, double xp, double yp, double delrx)
{
    double l2=pow(dist_bc(bc, hx[1]-hx[0], hy[1]-hy[0], fov[0], fov[1], delrx),2);
    if (l2==0) {
        return dist_bc(bc, hx[0]-xp, hy[0]-yp, fov[0], fov[1], delrx);
    }
    
    double tp=dot_bc(bc, xp-hx[0], yp-hy[0], hx[1]-hx[0], hy[1]-hy[0], fov[0], fov[1], delrx)/l2;
    
    if (tp<0) {
        return dist_bc(bc, hx[0]-xp, hy[0]-yp, fov[0], fov[1], delrx);
    }
    else if(tp>1.0){
        return dist_bc(bc, hx[1]-xp, hy[1]-yp, fov[0], fov[1], delrx);
    }
    else{
        double px=hx[0]+tp*(hx[1]-hx[0]);
        double py=hy[0]+tp*(hy[1]-hy[0]);
        return dist_bc(bc, px-xp, py-yp, fov[0], fov[1], delrx);
    }
}

//closest point on the link to point (xp, yp)
array<double,2> Link::get_intpoint(double xp, double yp)
{
    array<double,2> coordinates;
    double l2 = pow(hypot( hx[1]-hx[0], hy[1]- hy[0]) , 2);
    if (l2==0) {
        coordinates[0]=hx[0];
        coordinates[1]=hy[1];
    }
    double tp=dot(xp-hx[0],yp-hy[0],hx[1]-hx[0],hy[1]-hy[0])/l2;
    if (tp<0) {
        coordinates[0]=hx[0];
        coordinates[1]=hy[0];
    }
    else if(tp>1.0){
        coordinates[0]=hx[1];
        coordinates[1]=hy[1];
    }
    else{
        coordinates[0]=hx[0]+tp*(hx[1]-hx[0]);
        coordinates[1]=hy[0]+tp*(hy[1]-hy[0]);
    }
    return coordinates;
}

double Link::get_int_angle(double xp, double yp)
{
    double angle;
    double xcor,ycor;
    double slope=(hy[1]-hy[0])/(hx[1]-hx[0]);
    double yintercept = ycm - slope * xcm;
    xcor=(slope*yp + xp - slope*yintercept)/(slope*slope + 1);
    ycor=(slope*slope*yp + slope*xp + yintercept)/(1 + slope*slope);
    angle=atan2((ycor-yp),(xcor-xp));
    return angle;
}

vector<array<int, 2> > Link::get_quadrants()
{
    return quad;
}

array<double, 2> Link::get_direction()
{
    return {cos(phi), sin(phi)};
}
