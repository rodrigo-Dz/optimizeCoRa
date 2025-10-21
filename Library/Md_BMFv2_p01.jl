# Brink motif feedback (v01)
#   with Y inducing B synthesis,

# Julia v.1.8

module mm
	# Required libraries
	using DifferentialEquations
	using ParameterizedFunctions

    # ODE system with feedback
    function FB(du, u, p, t)
        Y, A, B, C, U, Us = u
        g, mY, gY, mA, kD, mB, mU, eP, e0, bA, bI, mAs = p
    
		du[Y]  = (mY * U)              - ((g + gY) * Y)
		du[A]  = (mA * (kD/(kD + Y)))  - (g * A) - (eP * A * B) + (e0 * C) - (bA * A * Us)
		du[B]  =  mB                   - (g * B) - (eP * A * B) + (e0 * C) - (bI * B * U)
		du[C]  =                       - (g * C) + (eP * A * B) - (e0 * C) + (bA * A * Us)
		du[U]  =  mU                   - (g * U) + (bA * A * Us) - (bI * B * U)
		du[Us] =                       - (g * Us) - (bA * A * Us) + (bI * B * U)
    end

    # ODE system without feedback
    function nFB(du, u, p, t)
        Y, A, B, C, U, Us = u
        g, mY, gY, mA, kD, mB, mU, eP, e0, bA, bI, mAs = p
    
		du[Y]  = (mY * U)              - ((g + gY) * Y)
		du[A]  =  mAs                  - (g * A) - (eP * A * B) + (e0 * C) - (bA * A * Us)
		du[B]  =  mB                   - (g * B) - (eP * A * B) + (e0 * C) - (bI * B * U)
		du[C]  =                       - (g * C) + (eP * A * B) - (e0 * C) + (bA * A * Us)
		du[U]  =  mU                   - (g * U) + (bA * A * Us) - (bI * B * U)
		du[Us] =                       - (g * Us) - (bA * A * Us) + (bI * B * U)
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
		p[:mBs] = p[:mB] * ss[1];
	end;
end