/***
* Name: arrays
* Author: allanlao
* Description: 
* Tags: Tag1, Tag2, TagN
***/
model arrays

global {
/** Insert the global definitions, variables and actions here */
	init {
		
		ask cell{
			grid_value <- 8;
		}
		
		do build_walls;
	}

	action build_walls {
		ask cell {

		//top
			if (grid_value >= 8.0) {
			//	walls << line(shape.points where (each.y < location.y));
				neighbours <- neighbours where (each.location.y >= (location.y));
			}
			//left
			if (not even(grid_value)) {
			//	walls << line(shape.points where (each.x < location.x));
				neighbours <- neighbours where (each.location.x >= (location.x));
			}
			//bottom
			if (grid_value = 2.0 or grid_value = 3.0 or grid_value = 6.0 or grid_value = 7.0 or grid_value = 10.0 or grid_value = 11.0 or grid_value = 14.0 or grid_value = 15.0) {
			//	walls << line(shape.points where (each.y > location.y));
				neighbours <- neighbours where (each.location.y <= (location.y));
			}
			//right
			if ((grid_value >= 4.0 and grid_value <= 7.0) or (grid_value >= 12.0 and grid_value <= 15.0)) {
			//	walls << line(shape.points where (each.x > location.x));
				neighbours <- neighbours where (each.location.x <= (location.x));
			}

	

			add self to: neighbours at: 0;
		}
				ask cell {
				loop c over: copy(neighbours) {
					if not (self in c.neighbours) {
						using topology(world) {
							neighbours <- neighbours where ((each.location distance_to c.location) >= (each.location distance_to location));
						}

					}

				}

			}

	}

}

grid cell neighbors: 8 width: 112 height: 27 {
	list<cell> neighbours;
}

experiment arrays type: gui {
/** Insert here the definition of the input and output of the model */
	output {
		display main{
			grid cell lines:#black;
		}
	}

}
