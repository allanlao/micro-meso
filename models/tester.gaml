/**
* Name: tester
* Author: allanlao
* Description: 
* Tags: Tag1, Tag2, TagN
*/

model tester

global {
	
	init{
	list<people> pips <- nil;
	write ""+  pips count (each.name);
	
	}
	
}


species people control:fsm{
	
	state born initial:true{
		write "in initial state";
	}
	state meso{
		write "in meso state";
		
	}
	state micro{
		write "in micro state";
	}
	
	reflex dothis{
		write "doing reflex";
	}
	
}


species meso_block{
	
}

species mesocell{
	list<cell> member_cell;
	point center;
	int density;
	
	
	aspect base{
		draw circle(20) at: center color:#red;
	}
	

}

grid cell neighbors:8 height:2 width:2{
	init{
		grid_value <- 1.0;
	}
	aspect base{
	 draw circle(1) ;	
	}
}

experiment tester type: gui {
	/** Insert here the definition of the input and output of the model */
	output {
		display main{
			grid cell lines:#black;
			species cell aspect:base;
			species mesocell aspect:base;
		}
	}
}
