/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set syntax=cpp fenc=utf-8 ft=tcl et sw=4 ts=4 sts=4

   This file is part of the Feel library

   Author(s): Vincent HUBER <vincent.huber@cemosis.fr>

   Date 2016-01-22

   Copyright (C) 2013 Université de Strasbourg

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
   */

#include <feel/feel.hpp>
using namespace Feel;
/*
 * dt u - a lap(u) + grad(p) = g
 *                    div(u) = 0
 */
int main(int argc, char**argv )
{
  po::options_description app_options( "MyBackend options" );
  app_options.add_options()
    ( "movePhi", Feel::po::value<bool>()->default_value( false ), "Move Phi ?" );

  app_options.add(backend_options("phi"));
  Environment env(_argc=argc, _argv=argv,
      _desc=app_options,
      _about=about(_name="nonSteadyStokes",
        _author="Vincent HUBER",
        _email="vincent.huber@cemosis.fr"));

  auto mesh = loadMesh(_mesh=new Mesh<Simplex< 2 > > );

  auto Vh = THch<1>( mesh );

  auto Xh = Pch<3>( mesh );

  auto U  = Vh->element("U");
  auto V  = Vh->element("V");
  auto u = U.element<0>("u");
  auto v = V.element<0>("v");
  auto p = U.element<1>("p");
  auto q = V.element<1>("q");

  auto phi   = Xh->element("phi");
  auto phi_m = Xh->element("phi");
  phi = vf::project(
      _space=phi.functionSpace(),
      _range=elements(mesh),
      _expr=expr(soption("functions.f"))
      );

  auto e = exporter( _mesh=mesh );

  auto M_bdf = bdf(_space = Vh);
  M_bdf->start();
  M_bdf->initialize(U);

  auto M_bdf_phi = bdf(_space = Xh,_name="phi");
  M_bdf_phi->start();
  M_bdf_phi->initialize(phi);

  auto l  = form1(Vh);
  auto ls = form1(Vh);
  ls = integrate(_range=elements(mesh),_expr= inner(expr<2,1>(soption("functions.g")),id(v)));

  auto as = form2( _trial=Vh, _test=Vh );
  auto a  = form2( _trial=Vh, _test=Vh );

  auto lp  = form1(Xh);
  auto ap  = form2(_trial=Xh, _test=Xh);


  // Static part

  do{
    auto M_bdf_poly = M_bdf->polyDeriv();
    l = ls;
    l += integrate(_range=elements(mesh),_expr = inner(idv(M_bdf_poly.element<0>()),id(v)));

    a = integrate(_range=elements(mesh),
        _expr = (doption("parameters.m")*chi(idv(phi)>0.)+ doption("parameters.n")*chi(idv(phi)<=0.))* inner(gradt(u), grad(v))
        -div(v)*idt(p)+divt(u)*id(q)
        ); 
    a+= integrate(_range=elements(mesh),
        _expr=inner(idt(u),id(v))*M_bdf->polyDerivCoefficient(0)
        );

    a+=on(_range=markedfaces(mesh,"null"), _rhs=l, _element=u,_expr=vec(cst(0.),cst(0.)));
    a+=on(_range=markedfaces(mesh,"g"), _rhs=l, _element=u,_expr=vec(cst(0.),cst(-1)));
    a+=on(_range=markedpoints(mesh,"p"), _rhs=l, _element=p, _expr=cst(0.)); // To avoid indefinite pressure

    a.solve(_rhs=l,_solution=U); // Compute with default backend
    M_bdf->shiftRight(U);
    M_bdf->next();
    /*
     * Update phi
     */
    if(boption("movePhi"))
    {
      LOG(INFO) << "Moving Phi\n";
      ap = integrate(_range=elements(mesh),
          _expr = M_bdf_phi->polyDerivCoefficient(0)*id(phi)*idt(phi)
          + inner(idv(u), trans(grad(phi)))*idt(phi));
      LOG(INFO) << "LHS Done\n";

      lp = integrate(_range=elements(mesh),
          _expr = idv(M_bdf_phi->polyDeriv())*id(phi));
      LOG(INFO) << "RHS Done\n";

      ap.solve(_rhs=lp,
          _solution=phi,
          _name="phi");

      LOG(INFO) << "Solve Done\n";

      M_bdf_phi->shiftRight(phi);
      M_bdf_phi->next();
    }

    e->step(M_bdf->time()) -> add( "u", u );
    e->step(M_bdf->time()) -> add( "p", p );
    e->step(M_bdf->time()) -> add( "phi", phi );
    e->save();
  }while(M_bdf->isFinished() == false);

}

