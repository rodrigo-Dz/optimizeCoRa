# Brink motif feedback (v01)
#   with Y inducing B synthesis,

# Julia v.1.8

module mm
	# Required libraries
	using DifferentialEquations
	using ParameterizedFunctions

    # ODE system with feedback
    function FB(du, u, p, t)
        Y, U, W, C = u
        g, mY, gY, mU, gU, mW,kD, gW, e0, eP, eM, mUs = p
    
		du[Y] = (mY * W)            - ((g + gY) * Y)
		du[U] = (mU * Y)            - ((g + gU) * U) - (eP * U * W) + ((e0 + gW + eM) * C)
		du[W] = (mW * (W/(W + kD))) - ((g + gW) * W) - (eP * U * W) + ((e0 + gU) * C)
		du[C] =                     - ((g + eM) * C) + (eP * U * W) - ((gU + gW + e0) * C)
    end

    # ODE system without feedback
    function nFB(du, u, p, t)
        Y, U, W, C = u
        g, mY, gY, mU, gU, mW,kD, gW, e0, eP, eM, mUs = p
    
		du[Y]  = (mY * W)            - ((g + gY) * Y)
		du[U]  =    mUs              - ((g + gU) * U) - (eP * U * W) + ((e0 + gW + eM) * C)
		du[W]  = (mW * (W/(W + kD))) - ((g + gW) * W) - (eP * U * W) + ((e0 + gU) * C)
		du[C]  =                     - ((g + eM) * C) + (eP * U * W) - ((gU + gW + e0) * C)
    end
    

	# Define system's output (total Y):
	function outFB(ss)
		return ss[1];
	end;
	function outNF(ss)
		return ss[1];
	end;

	# Define locally analogous system:
	function localNF(p,ss)
		p[:mUs] = p[:mU] * ss[1];
	end;
end