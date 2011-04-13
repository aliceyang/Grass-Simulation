Alice Yang
CIS 563 Homework 1
Jello Simulation
README

Questions
1.	Ks is the spring constant-- it determines how stiff the jello is (the higher the
	value of Ks, the stiffer the jello). Kd is the dampening constant. Kd is the 
	damping constant. It reduces the amount of oscillation on the jello.
	
2.	One of the benefits of this collision system is that the user can choose the integration
	type to use. One drawback is that the spring constants are limited by the size of the
	time steps in the integration methods. If the constnats are too big, the system blows up.
	
3.	One example of a stiff constrained system is a bead on a wire.

4.	While explicit integration methods advances the system forward in time, implicit integration
	evaluates the ODE at the point we're aiming at rather than where we came from. In other words,
	it advances the system backwards in time.
	
5.	You shouldn't use a penalty constraint to hold a bead on a wire because weak springs would
	lead to goopy constraints, while strong strings would propel the bead off its original
	trajectory. Instead, you should use differential constraints (implicit)
	
6.	The jello behaved realistically to an extent (bouncing vertically off the groud looked fairly real). 
	However, in more complicated cases the behvaior was off (ie. when bouncing off a slanted cylinder
	the bottom jiggled too much).
	
Design Decisions:

	I found out that whether or not the jello exploded depended largely on the constants I chose. 
	Therefore, even though some of the constants may look off, they were good enough to keep everything
	from exploding. 
	
	One problem I had was the cylinder endcap behavior. When the jello is too close to an endcap, sometimes
	it will exhibit erratic behavior. I've tried to rectify this by checking if t is close to 0 or 1 at the
	point of contact for each particle, and setting the normal and distance accordingly, but it doesn't
	fix all cases. Also, I decreased the time steps for Euler and Midpoint integration to stabilize the
	system.
	
	I discussed ideas about the assignment with Nop and Marley.
	I also found some slides online that were helpful in answering the questions: www.cs.cmu.edu/~15869-f10/lec/04/lec04.pdf