#include <regs.h>
#include <Eigen/LU>
#include <mp_nonlinear_odes.h>


mp_mat mpFindA(const mp_vec& t, const mp_mat& U, int m)
{
	mp_mat A(m,m);
	
	for(int i = 0; i<m; i++){
		for(int j = 0; j<m; j++){
				A(i,j) = mpInnerProd(U.row(i), U.row(j), t);
		}
	}
	return A;
}

mpreal mpPsuedoDet(const mp_mat& A)
{
	mpreal g = "1";
	int m = A.cols();
	for(int i=0; i<m; i++){
		g *= A(i,i);
	}
	return g;
}

mp_vec mpFindP(const mp_vec& t, const mp_mat& U, const mp_vec& dx, int m)
{
	mp_vec P(m);
	
	for(int i = 0; i<m; i++)
	{
		P(i) = mpInnerProd(U.row(i), dx, t);
	}
	
	return P;
}

mpreal mpInnerProd(const mp_vec& u1, const mp_vec& u2, const mp_vec& time)
{
	int mid = time.size();
	int n = u1.size()/mid;
	
	mp_vec a = u1.array() * u2.array();
	
	mp_vec aij = a.head(mid);
	
	if(mid != n*mid){
		for(int k=1; k<n; k++){
			aij += a.segment(k*mid, mid);
		}
	}
	
	return mpSimpson(time, aij);
}

vec regularization(soln_env *env){
	mpreal::set_default_prec(BITS);
	
	int n = (*env->measurements).rows(); 
	int m = (*env->initial_params).size(); 
	int lt = (*env->time).size(); 		
	
	mp_vec uNot(m);
	uNot = (*env->initial_params).cast<mpreal>();
	mp_mat bob(n*(m+1), lt);
	mp_mat U(n, m*lt);
	mp_mat A(m,m);
	mp_mat AT(m,m);
	mp_mat B = mp_mat::Identity(m, m);
	/* for(int i=0; i<m-1; i++){
		B(i+1,i) = -1;
	} 
	B = B.transpose()*B;*/
	mp_vec P(m);
	mp_vec du(m);
	mp_vec du1(m);
	du = uNot.cast<mpreal>();
	mpreal gamma;
	
	mp_vec t(lt); t << (mp_vec)(*env->time).cast<mpreal>();
	mp_mat soln(n,lt); soln << *env->mp_nth_soln;
	mp_mat msmt(n,lt); msmt << *env->mp_measurements;
	
	int LIMIT = 250;
	for(int i = 0; i<LIMIT; i++)
	{
		bob = mpQLinearRungeKutta4(*env->ode, t, uNot, (*env->initial_cond).cast<mpreal>(), soln);
		U = mpReshape(bob.bottomRows(n*m), m, n*lt);
		soln = bob.topRows(n);
		
		A = mpFindA(t, U, m); note(mp_cond(A));
		AT = A.transpose();
		P = mpFindP(t, U, mpReshape(msmt-soln, 1, n*lt).row(0), m);
		
		//du = (mpCofactor(A)/mpPsuedoDet(A))*P;
		
		gamma = "2";
		for(int j=0; j<1500; j++)
		{
			//du = gamma*du;
			//du = (AT*A + gamma*gamma*B).householderQr().solve(AT*P);
			du = mp_inverse(AT*A + gamma*gamma*B)*(AT*P); 
			//du = mp_inverse(A)*P; 
			
			gamma *= ".85";
			note(mp_norm(du));
			note(mp_norm(A*du-P));
			if(mp_norm(A*du-P) < 0.1 ){
				break;
			}
			if(mpfr::isnan(du.norm())){
				break;
			}
		}

		uNot += du;
		note(uNot);
		if(mp_norm(du) < 1E-4  || mpfr::isnan(du.norm())){
			 break;
		} else if (i >= LIMIT-1){
			log_err("Function did not converge.");
			note("u = ");
			note(uNot);
			exit(1);
		}
	}
	return uNot.cast<double>();
}

mp_mat mpQLinearRungeKutta4(string fname, const mp_vec& time, const mp_vec& u, const mp_vec& yNot, const mp_mat& xNminus)
{
	thesis::mpnonlinearOdes no;
	mp_sys fhandle = no.mpOdeFuncMap[fname];

	int N = time.size();
	int n = xNminus.col(1).size();
	
	//number of equations
	int m = yNot.size();
	
	//timestep
	mpreal h = time(1)-time(0);
	
	//init stuff
	mp_mat w(m, N);
	w.fill(0);
	w.col(0) = yNot;
	
	mp_vec k1(m), k2(m), k3(m), k4(m);
	
	thesis::mp_spline Xn[n];
	for(int i = 0; i < n; i++){
		Xn[i].update(time, xNminus.row(i));
	}
	
	for (int i = 0; i<N-1; i++)
	{
		
		k1 = h*mpQlinear(fhandle, time(i),       w.col(i),        u, Xn, n);
		k2 = h*mpQlinear(fhandle, time(i) + h/2, w.col(i) + k1/2, u, Xn, n);
		k3 = h*mpQlinear(fhandle, time(i) + h/2, w.col(i) + k2/2, u, Xn, n);
		k4 = h*mpQlinear(fhandle, time(i) + h,   w.col(i) + k3,   u, Xn, n);
		
		//cout <<"(" << i <<")\nk1"<< k1 << "\nk2:" << k2 << "\nk3:" << k3 << "\nk4:" << k4 <<endl;
		w.col(i+1) = w.col(i) + (k1 + 2*(k2 + k3) + k4)/6;
	}
	
	return w;
}

mpreal mpSimpson(const mp_vec& t, const mp_vec& x)
{
	thesis::mp_spline sim(t,x);
	mpreal temp;
	
	int N = 1501;
	int end = t.size()-1;
	mpreal h = (t(end)- t(0))/(N-1);
	
    mpreal area = x(0) + x(end);
    for(int i = 1; i<N-1; i++)
	{
		temp = sim.interpolate(t(0)+h*i);
        if((i+1)%2 == 0)
		{
            area += 2*temp;
        }
		else
		{
            area += 4*temp;
		}
    }
	
    area = h/3 * area;
	return area;
}

mp_mat mpDer(const mp_mat& dx, const mpreal& dt){
	mp_mat ans(dx.size(),1);
	ans << dx/dt;
	return ans;
}

mp_mat mpJac(mp_sys f, mpreal t, const mp_mat& x, const mp_mat& u, const mpreal& h){
	int n = x.size();
	mp_mat dx(n,n);
	dx << mp_mat::Identity(n,n)*h;
	mp_mat fprime(n,n);
	for(int j=0; j<n; j++)
	{
		fprime.col(j) = -f(t, x+2*dx.col(j), u) + 8*f(t, x+dx.col(j), u) - 8*f(t, x-dx.col(j), u) + f(t, x-2*dx.col(j), u);
	}
	
	return fprime/(12*h);
}

mp_mat mpQlinear(mp_sys fhandle, const mpreal& t, const mp_vec& x, const mp_vec& u, thesis::mp_spline* Xn, int n)
{
	int m = u.size();
	mpreal step = "2.2E-16";
	
	mp_vec xn1(n); // this is x_N-1
	mp_vec dxn(n);
	
	for(int i=0; i<n; i++)
	{
		xn1(i) = Xn[i].interpolate(t);
		dxn(i) = x(i) - xn1(i); // x_N - x_{N-1}
	}
	
	mp_mat dx(n,n);// for x derivative
	dx << mp_mat::Identity(n,n)*step;
	
	//First few lines of linearization
	mp_mat fx(n,1);
	fx = fhandle(t, xn1, u);
	mp_mat ans(n+n*m, 1);
	ans << mp_mat::Zero(n+n*m, 1);
	ans.block(0, 0, n, 1) = fx + mpJac(fhandle, t, xn1, u, step)*dxn;

	/* 	
	ans.block(0, 0, n, 1) = fx;
	for(int j=0; j<n; j++)
	{
		ans.block(0, 0, n, 1) += mpDer(fhandle(t, xn1+dx.col(j), u) - fx, step)*dxn(j);
	} 
	*/
	
	mp_mat dfdx(n,1);
	mp_mat dfdu(n,1);
	
	// the Un part of t the linearization
	mp_mat dun(m, m);// for u derivative
	dun << mp_mat::Identity(m,m)*step;
	
	int ind;
	for(int j=0; j<m; j++)
	{
		ind = (j+1)*n;
		dfdu = fhandle(t, xn1, u+dun.col(j));
		ans.block(ind, 0, n, 1) = mpDer(dfdu - fx, step); // df/du
		for(int k=0; k<n; k++){
			dfdx = fhandle(t, xn1+dx.col(k), u);
			ans.block(ind, 0, n, 1) += mpDer(dfdx  - fx, step)*x(ind+k); //J*Un
			//ans.block(ind, 0, n, 1) += mpDer(fhandle(t, xn1+dx.col(k), u+dun.col(j)) - dfdu - dfdx + fx, step*step)*dxn(k); //phi_ij
			ans.block(ind, 0, n, 1) += mpDer(
				  fhandle(t, xn1+dx.col(k), u+dun.col(j)) 
				- fhandle(t, xn1+dx.col(k), u-dun.col(j)) 
				- fhandle(t, xn1-dx.col(k), u+dun.col(j))
				+ fhandle(t, xn1-dx.col(k), u-dun.col(j)), 4*step*step)*dxn(k); //phi_ij
		}
	}
	
	return ans;
}

mpreal mp_norm(const mp_mat& M){
	return M.norm();
}

mp_mat mp_inverse(const mp_mat& M){
	return M.inverse();
}

mpreal mp_cond(const mp_mat& A){
	return A.norm()*A.inverse().norm();
}



mp_mat mpReshape(const mp_mat& U, int n, int m)
{
	mp_mat newU(n,m);
	newU.fill(0);
	int olt = U.row(0).size(); 	//old time(row) length
	int on = m/olt;				//on col length
	
	for(int i = 0; i<n; i++){
		for(int j = 0; j<on; j++){
			newU.block(i, j*olt, 1, olt) = U.row(i*on + j);
		}
	}
	return newU;
	
	
}

mp_mat mpCofactor(const mp_mat& A){
	//check it is square
	int m = A.cols();
	mp_mat C;
	C = A;
	for(int i=0; i<m; i++){
		for(int j=0; j<m; j++){
			C(i,j) = pow(-1,i+j)*mpRowColRemoval(A,i,j).determinant();
		}
	}
	return C;
} 

mp_mat mpRowColRemoval(const mp_mat& A, const int row, const int col){
	int mm = A.cols()-1;
	mp_mat P(mm,mm);
	
	P.topLeftCorner(row, col) = A.topLeftCorner(row, col);     // P(1:rows, 1:cols)
	P.topRightCorner(row, mm-col) =	A.topRightCorner(row, mm-col);    // P(1:rows, end-cols+1:end)
	P.bottomLeftCorner(mm-row, col)  = A.bottomLeftCorner(mm-row, col);  // P(end-rows+1:end, 1:cols)
	P.bottomRightCorner(mm-row, mm-col) = A.bottomRightCorner(mm-row, mm-col); // P(end-rows+1:end, end-cols+1:end)
	
	return P;
}/* */