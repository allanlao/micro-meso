model micromeso
//test branch stay 
import "../includes/common_pedestrian.gaml"


global
{
	string folder <- 'new_bldg';
	//	string folder <- 'nb101';
	file micro_layer <- grid_file("../includes/paper/" + folder + "/micro_layer.asc");
	file meso_layer <- grid_file("../includes/paper/" + folder + "/meso_layer.asc");
	file basin_file <- csv_file("../includes/paper/" + folder + "/basins.csv",",");
	file trap_file <- file("../includes/paper/" + folder + "/trap.csv");
	file swamp_file <- csv_file("../includes/paper/" + folder + "/swamp.csv");
	file counter_file <- file("../includes/paper/" + folder + "/counters.csv");
	//list<cell> basin_cells <- nil;
	list<int> basin_group <- nil;
	list<cell> free_cell <- shuffle(cell);
	//list<zombie> zombies <- nil;
	float egress_time <- 0.0;
	int total_population <- 75;
	string door <- "right";
	string smethod <- "tabu";
	map egress <- nil;
	bool load_navigation <- false;
	int simulation_run <- 3;
	float exp_time_start <- 0.0;
	float exp_time_end <- 0.0;
	int people_ctr <- 0;
	list<int> counters_list <- nil;
	geometry shape <- envelope(micro_layer);
	bool showGrid <- true;
	list<cell> source_cell <- nil;
	map prob_dist_map <- nil;
	int dist_k <- 200;
	int dist_n <- 0;
	float dist_p <- 0.5;
	int flow <- 3;
	int W <- 0;
	int L <- 0;
	int pop_counter <- 10;
	bool save_ntxy <- false;
	bool save_zombies <- false;
	bool showResult <- false;
	bool save_result <- true;
	float actual_egress_time <- 24.0;
	list<list<float>> fundamental_list <- nil;
	bool stay_in_place <- true;
	float max_spv <- 1.0;
	//Oranges, Greens, Blues, Greys
	// list<rgb> colr <- brewer_colors("Greens", max_density); 
	 list<rgb> colr <- ([rgb(247,247,247,255),rgb(204,204,204,255),rgb(150,150,150,255),rgb(99,99,99,255),rgb(37,37,37,255)]);
	 
    int penalty <- 0;
	int default_reward <- 0;
	string system_result_filename <- "/results/system/results.csv";
	
	float micro_cell_diameter <- 0.5;
	
	
	
	list<int> meso_group <- remove_duplicates(list<int>(meso_cell collect (each.grid_value)));

	init
	{
		write "initalizing : " + total_population;
	
write basin_file;
		// "dont forget the trap so waypoints will work";
		//max_density <- int((cell_diameter >= 1) ? ceil(((cell_diameter / 0.5) ^ 2) + 1) : ceil((cell_diameter / 0.5) ^ 2));
		do build_trap;
		do build_walls;
		do build_basins;
		do build_rcso_counters;
		do build_swamps;
		if (load_navigation)
		{
			do loadNavigation;
		} else
		{
			do build_navigation;
		}
       
      
        
       
       	loop g over:meso_group{
       		write g;
       		if ( g >0){
				create meso_block{
					member_cells <- meso_cell where(each.grid_value = g);
					write member_cells;
					center <- point((member_cells[0].location.x + member_cells[3].location.x)/2, (member_cells[0].location.y + member_cells[3].location.y)/2);	
				}
			}
       }
       
         ask cell{
         	if (meso_cell[grid_x,grid_y].grid_value > 0){
         		is_meso <- true;
         	}
        	do update_color;
        	
        }  
       
      ask spvcells{
      
      	walls <- cell[grid_x,grid_y].walls;
       	spv <- cell[grid_x,grid_y].spv;
       	do update_color;
      }
       
        //list of micro cell that needs to be updated
      
		list<cell> cell_no_basin <- cell where (each.basin = 0);
		source_cell <- cell where (each.basin = 2);
		source_cell <- source_cell sort (each.basin);

		//add distribution to source cells;
		dist_n <- length(source_cell) - 1;
		prob_dist_map <- basin_map_probability(dist_n, dist_k, dist_p);
		
		if (!use_people_generator)
		{
			source_cell <- (cell - (cell where (each.basin = 1))) where (each.is_trap);
			
			//source cells can only be micro cells
			create people number: total_population
			{
				my_cell <- one_of(source_cell where (length(people inside (each)) =0));
				location <- my_cell.location;
				nav_id <- 1; //should be from destination list
				time_start <- float(cycle + 1);
				if (my_cell.is_meso){
				position_in_cell <- (cell_diameter / 2);
				//if people is created at once then initial position should change ;
				locations_in_cell << [(cycle) * delta_time, position_in_cell];
				entry_angle <- 0;
				just_in <- false;
				}else{
				position_in_cell <- (cell_diameter / 4);
				//if people is created at once then initial position should change ;
				locations_in_cell << [(cycle) * delta_time, position_in_cell];
				entry_angle <- 0;
				just_in <- false;
				}
				
			}

		}


	 	
		do update_density;
		ask updatable_micro_space{
			write "drawing circles";
			draw circle(100)  color:#red;
		} 
		
	} //end of init
	
	
     list<cell> updatable_micro_space <- cell where (!each.is_meso and  (each.neighbours != nil));
	
	//my functions
	reflex create_people when: use_people_generator and (people_ctr < total_population)
	{
	//run generator
	
	loop times:flow{
		
	
		loop x over: prob_dist_map.keys
		{
			x <- int(x);
			if (length(people inside source_cell[x]) < max_density)
			{
				bool head <- flip((float(prob_dist_map[x]) / dist_k * 50 + 50) / 100);
				if (head and people_ctr < total_population)
				{
					create people number: 1
					{
						my_cell <- source_cell[x];
						location <- my_cell.location;
						nav_id <- 1; //should be from destination list
						time_start <- float(cycle + 1);
						
						if (my_cell.is_meso){
							position_in_cell <- (max_speed * delta_time) with_precision 2;
							locations_in_cell << [cycle * delta_time, 0];
							locations_in_cell << [(cycle + 1) * delta_time, position_in_cell];
						
						entry_angle <- source_cell[x].entry_angle;
						}else{
							position_in_cell <- (cell_diameter /4) with_precision 2;
							locations_in_cell << [cycle * delta_time, 0];
							locations_in_cell << [(cycle + 1) * delta_time, position_in_cell];
						
						entry_angle <- source_cell[x].entry_angle;
						}
					}

					people_ctr <- people_ctr + 1;
				}

			}

		}
		
		}

		do update_density;
		
	}

  
    


	//my actions
	action build_walls
	{
		ask cell
		{

		//top
			if (grid_value >= 8.0)
			{
				walls << line(shape.points where (each.y < location.y));
				neighbours <- neighbours where (each.location.y >= (location.y));
			}
			//left
			if (not even(grid_value))
			{
				walls << line(shape.points where (each.x < location.x));
				neighbours <- neighbours where (each.location.x >= (location.x));
			}
			//bottom
			if (grid_value = 2.0 or grid_value = 3.0 or grid_value = 6.0 or grid_value = 7.0 or grid_value = 10.0 or grid_value = 11.0 or grid_value = 14.0 or grid_value = 15.0)
			{
				walls << line(shape.points where (each.y > location.y));
				neighbours <- neighbours where (each.location.y <= (location.y));
			}
			//right
			if ((grid_value >= 4.0 and grid_value <= 7.0) or (grid_value >= 12.0 and grid_value <= 15.0))
			{
				walls << line(shape.points where (each.x > location.x));
				neighbours <- neighbours where (each.location.x <= (location.x));
			}

			ask cell
			{
				loop c over: copy(neighbours)
				{
					if not (self in c.neighbours)
					{
						using topology(world)
						{
							neighbours <- neighbours where ((each.location distance_to c.location) >= (each.location distance_to location));
						}

					}

				}

			}

			add self to: neighbours at: 0;
		}

	}

	action build_trap
	{
		int x1 <- 0;
		int x2 <- 0;
		int y1 <- 0;
		int y2 <- 0;
		matrix m <- matrix(trap_file);
		x1 <- int(m[0, 0]);
		y1 <- int(m[1, 0]);
		x2 <- int(m[2, 0]);
		y2 <- int(m[3, 0]);
		loop c from: x1 to: x2
		{
			loop r from: y1 to: y2
			{
				cell[c, r].is_trap <- true;
			}

		}

	}
	
	
	

	action build_basins
	{
	//col,row,group,type,reward
		matrix m <- matrix(basin_file);
		
		write sample(m);
		
		//skip header
		loop row_ctr from: 0 to: m.rows - 1
		{
			

			int x1 <- 0;
			int x2 <- 0;
			int y1 <- 0;
			int y2 <- 0;
			x1 <- int(m[0, row_ctr]);
			y1 <- int(m[1, row_ctr]);
			x2 <- int(m[2, row_ctr]);
			y2 <- int(m[3, row_ctr]);
			loop col from: x1 to: x2
			{
				loop row from: y1 to: y2
				{
					cell[col, row].basin_group <- int(m[4, row_ctr]);
					cell[col, row].basin <- int(m[5, row_ctr]);
					cell[col, row].reward <- int(m[6, row_ctr]);
					
				}

			}

		}

	}
	
	
	
	action build_swamps
	{
	//col,row,group,type,reward
	 
		matrix m <- matrix(swamp_file);
		
		write "building swamp";
		write m.rows;
	
		if (m.rows <= 1)
		{
			
			return;
		}
		
		loop row_ctr from: 0 to: m.rows - 1
		{
			if (m[0, row_ctr] = "#")
			{
				return;
			}

			int x1 <- 0;
			int x2 <- 0;
			int y1 <- 0;
			int y2 <- 0;
			x1 <- int(m[0, row_ctr]);
			y1 <- int(m[1, row_ctr]);
			x2 <- int(m[2, row_ctr]);
			y2 <- int(m[3, row_ctr]);
			loop col from: x1 to: x2
			{
				loop row from: y1 to: y2
				{
					cell[col, row].is_swamp <-true;
					
				}

			}

		}

	}

	action build_rcso_counters
	{
	//col,row,group,type,reward
		matrix m <- matrix(counter_file);
		if (m.rows <= 1)
		{
			return;
		}

		loop row_ctr from: 0 to: m.rows - 1
		{
			int x1 <- 0;
			int x2 <- 0;
			int y1 <- 0;
			int y2 <- 0;
			x1 <- int(m[0, row_ctr]);
			y1 <- int(m[1, row_ctr]);
			x2 <- int(m[2, row_ctr]);
			y2 <- int(m[3, row_ctr]);
			loop col from: x1 to: x2
			{
				loop row from: y1 to: y2
				{
					cell[col, row].rcso_counter_no <- int(m[4, row_ctr]);
					cell[col, row].color <- # orange;
					cell[col, row].is_counter <- true;
				}

			}

			add int(m[4, row_ctr]) to: counters_list;
		}

	}

	action loadNavigation
	{
		write "loading navigation";
		basin_group <- remove_duplicates(cell collect each.basin_group);
		basin_group <- basin_group where (each > 0);
		loop bg over: basin_group
		{
			file nav_file <- csv_file("../includes/paper/" + folder + "/navigation/" + bg + "_nav.csv", ",");

			
			matrix qm <- matrix(nav_file);
			
			//max_spv <- max(qm);

			//center
//			loop row from: 0 to: qm.rows - 1
//			{
//				loop col from: 0 to: qm.columns - 1
//				{
//				//get center
//					
//					cell[col, row].spv <- qm[col, row];
//				}
//
//			}
			add (bg)::qm to: nav_map;
		}

	}

	action build_navigation
	{
		
	//sink cells = 1
	//basin_cells <- cell where (each.basin = 1);
		basin_group <- remove_duplicates(cell collect each.basin_group);
		basin_group <- basin_group where (each > 0);
		loop bg over: basin_group
		{
			int iter <- 0;
			int cols <- cell max_of (each.grid_x) + 1;
			int rows <- cell max_of (each.grid_y) + 1;
			list<list<float>> q_list <- [];
			//	list<cell> von_neighbours <- [];
			list<list<cell>> state_action_list <- [];
			matrix reward_matrix <- default_reward as_matrix ({ cols, rows });
			matrix q_matrix <- 0.0 as_matrix ({ cols, rows });

			//cells of the same group
			//	list<cell> basins <- cell where ((int(each.basin / 10)) = int(od_matrix[1, od_ctr])); //od_matrix col1(destination) row0
			list<cell> basins <- cell where (each.basin_group = bg);
			float current_sum <- 0.0;
			float prev_sum <- 0.0;
			int counter <- 0;

			//build reward
			
			//assign penalty to walls
			
			loop g over:cell{
				int p <- 0;
					if( length(g.neighbours) < 9){
						p <- penalty;
					}
					
					if (g.is_swamp){
						p <- p + penalty;
					}
					
					reward_matrix[g.grid_x, g.grid_y] <- -p;
			}
			
		
			
			
			//get reward of basins from rewards file
			loop g over: basins
			{
				reward_matrix[g.grid_x, g.grid_y] <- g.reward;
			}

			//build qmatrix
			//von neuman
			ask cell
			{

			//	add neighbours collect (reward_matrix[each.grid_x, each.grid_y]) to: q_list1;
				add neighbours collect (0.0) to: q_list;
				add neighbours to: state_action_list;
			}

			//iterate up to max_iteration  or until diffence between previous and next state is nil
			loop i from: 0 to: max_iteration
			{

			//state_action_list is the neighborhood
				loop r from: 0 to: length(state_action_list) - 1
				{
					cell current_state <- state_action_list[r][0];
					//loop for von neuman
					loop c from: 0 to: length(state_action_list[r]) - 1
					{
					//q(s,a) per action and determine the maximum
						cell next_state <- state_action_list[r][c];
						float cd <- current_state.location distance_to next_state.location;
						float reward <- float(reward_matrix[next_state.grid_x, next_state.grid_y]);
					
						if (cd > 1.0)
						{
							reward <- sqrt(abs(reward));
						}

						float max_q <- max(q_list[int(next_state)] collect (each)); //what is the max reward
						q_list[r][c] <- floor(reward + learning_rate * max_q);
					}

					//then loop for moore - von neuman

				} //end of i state action cycle
				loop g over: q_list
				{
					current_sum <- sum(g collect (each));
				}

				if (abs(prev_sum - current_sum) <= 0.000001)
				{
					counter <- counter + 1;
				} else
				{
					counter <- 0;
				}

				if (counter > 100)
				{
					write "breaking qlearning iteration";
					break;
				}

				prev_sum <- current_sum;
				iter <- i;
			} //max_iter




			//convert q_list to navigation matrix
			int ctr <- 0;

			//center
			loop row from: 0 to: q_matrix.rows - 1
			{
				loop col from: 0 to: q_matrix.columns - 1
				{
				//get center
					q_matrix[col, row] <- q_list[ctr][0];
					ctr <- ctr + 1;
				}

			}

			//max_spv <- max(q_matrix);

			//center
			loop row from: 0 to: q_matrix.rows - 1
			{
				loop col from: 0 to: q_matrix.columns - 1
				{
				//get center
					q_matrix[col, row] <- q_matrix[col, row];
					cell[col, row].spv <- q_matrix[col, row];
				}

			}

		
			add bg::q_matrix to: nav_map;
			do save_matrix(q_matrix, "../includes/paper/" + folder + "/navigation/" + bg + "_nav.csv");
		} //build navs
	
	}

	action save_matrix (matrix m, string f)
	{
		save m to: f rewrite: true type: "csv" header: "false";
	}

	reflex remove_dead_people
	{
      ask people
		{
			if (is_dead and can_move_out() and my_cell.is_meso)
			{
				float speed <- S(length(people inside my_cell) - 1);
				position_in_cell <- (cell_diameter) with_precision 2;
				if (position_in_cell > cell_diameter / 2)
				{
				//check if midpoint exist if not add the midpoint in locations in cell
					float last_loc <- (last(locations_in_cell))[1];
					float last_time <- (last(locations_in_cell))[0];
					if (last_loc < cell_diameter / 2)
					{
						float new_d <- ((cell_diameter / 2) - last_loc) with_precision 2;
						locations_in_cell << [last_time + (new_d / speed) with_precision 1, cell_diameter / 2];
					}

				}

				if (my_cell.is_trap and my_cell.is_meso)
				{
					locations_in_cell << [(cycle + 1) * delta_time, position_in_cell];
					loop g over: locations_in_cell
					{
						float r <- 0.0;
						if (g[1] <= cell_diameter / 2)
						{
							r <- (cell_diameter / 2 - g[1]) with_precision 2;
							float x <- getXCoordinate(my_cell.grid_x, entry_angle, r);
							float y <- getYCoordinate(-my_cell.grid_y, entry_angle, r);
							waypoints << [g[0], x, y];
						} else
						{
							r <- (g[1] - cell_diameter / 2) with_precision 2;
							float x <- getXCoordinate(my_cell.grid_x, exit_angle, r);
							float y <- getYCoordinate(-my_cell.grid_y, exit_angle, r);
							waypoints << [g[0], x, y];
						}

					}

				}else{
					
				}

				create zombie
				{
				
			
					peoplename <- myself.name;
					waypoints <- myself.waypoints;
					time_expired <- last(waypoints)[0];
					cycles_expired <- cycle + 1;
					create performance_record number: 1 returns: created_performance_record;
					my_performance_record <- first(created_performance_record);
					rcso_counter_no <- myself.rcso_counter_no;
					exit_cell <- myself.my_cell.name;
				}
              
				do die;
			
				
			}

		}
		
		

	}

	reflex mark_people_destination
	{
		ask cell
		{
			temp_cell_density <- density;
		}

		ask people
		{
			do mark_next_destination;
		}

		//do pause;

	}

	reflex move_agents
	{
		ask people
		{
			do move;
		}

	}

	reflex update_micro_space
	{
		
		ask updatable_micro_space{
			draw circle(1)  color:#red;
		}
		
	}

	/* 	reflex update_fundamental_list{
		//list<people> <- people where (each.my_cell.is_trap);		
		
		int n <- people count (each.my_cell.is_trap);
		float d <- n/ (trap_area);
		// float speed <- S(my_cell.density - 1);
		
		fundamental_list << [(cycle * delta_time),d,0];
	
	//	write "total agent is " + n + "length >" + trap_length + " width " + trap_width + "density> " + d;
		
	} */
	reflex terminate when: empty(people)
	{
		do compute_performance;
		if (save_ntxy)
		{
			do save_ntxy;
		}

		if (save_zombies)
		{
			do save_zombies;
		}

		//do plot_egress;
		if (showResult)
		{
			do show_result;
		}

		write "terminating experiment";
		do die;
	}

	action show_result
	{
		write "Population:" + total_population;
		write "parameters: NavBias>" + Navigation_bias;
		write "egress :" + egress_time + " att:" + ave_travel_time + " atd:" + ave_walking_distance + " ave walking speed:" + ave_walking_speed;
	}

	action save_ntxy
	{
		string filename <- "/results/" + total_population + "_ntxyfile.csv";
		save ["n", "t", "x", "y", "p"] to: "../includes/paper/" + folder + filename type: "csv" header: false rewrite: true;
		ask zombie
		{
			loop ctr from: 0 to: length(waypoints) - 1
			{
				string line <- name + "," + waypoints[ctr][0] + "," + waypoints[ctr][1] + "," + waypoints[ctr][2] + "," + total_population;
				save (line) to: "../includes/paper/" + folder + filename type: "csv" header: false rewrite: false;
			}

		}

	}

	action compute_performance
	{
		ask zombie
		{
			my_performance_record.travelTime <- 0.0;
			//ind walking distance
			loop ctr from: 0 to: length(waypoints) - 2
			{

			// chebyshev
			// float d <- point(waypoints[ctr][1], waypoints[ctr][2]) distance_to point(waypoints[ctr + 1][1], waypoints[ctr + 1][2]);
				float d <- max([abs(waypoints[ctr][1] - waypoints[ctr + 1][1]), abs(waypoints[ctr][2] - waypoints[ctr + 1][2])]);

				//euclidean
				// float d <- point(waypoints[ctr][1], waypoints[ctr][2]) distance_to point(waypoints[ctr + 1][1], waypoints[ctr + 1][2]);
				my_performance_record.insWalkingDisplacement << d;
				my_performance_record.insWalkingDistance << abs(d);
				my_performance_record.insVelocity << (1 / delta_time) * d;
				my_performance_record.insSpeed << abs(d) / (abs(waypoints[ctr][0] - waypoints[ctr + 1][0]) > 0 ? abs(waypoints[ctr][0] - waypoints[ctr + 1][0]) : 0.01);
			}

			my_performance_record.travelTime <- (last(waypoints)[0] - first(waypoints)[0]) with_precision 2;
			my_performance_record.walkingDistance <- (sum(my_performance_record.insWalkingDistance)) with_precision 2;
			my_performance_record.speed <- (mean(my_performance_record.insSpeed)) with_precision 2;
			//my_performance_record.speed <- my_performance_record.walkingDistance / my_performance_record.travelTime;
			my_performance_record.walkingDisplacement <- (sum(my_performance_record.insWalkingDisplacement)) with_precision 2;
			my_performance_record.straightLineDistance <- 1 + (point(first(waypoints)[1], first(waypoints)[2]) distance_to point(last(waypoints)[1], last(waypoints)[2])) with_precision 2;
			my_performance_record.straightDistanceRatio <- (my_performance_record.walkingDistance / my_performance_record.straightLineDistance) with_precision 2;
			my_performance_record.velocity <- (sum(my_performance_record.insVelocity) / (length(my_performance_record.insVelocity)) / delta_time) with_precision 2;
			my_performance_record.walkingDelay <- ((my_performance_record.walkingDistance - my_performance_record.straightLineDistance) / my_performance_record.speed) with_precision 2;
			my_performance_record.straightLineSpeed <- (my_performance_record.straightLineDistance / ((length(waypoints) - 1) * delta_time)) with_precision 2;
			my_performance_record.movingDirection <- (my_performance_record.velocity / my_performance_record.straightLineSpeed) with_precision 2;
		}

		//get system averages
		ave_walking_distance <- (sum(zombie collect (each.my_performance_record.walkingDistance)) / length(zombie)) with_precision 2;
		ave_walking_speed <- (sum(zombie collect (each.my_performance_record.speed)) / length(zombie)) with_precision 2;
		ave_travel_time <- (sum(zombie collect (each.my_performance_record.travelTime)) / length(zombie)) with_precision 2;
		egress_time <- max(zombie collect (each.time_expired));
		int total_zombie <- length(zombie);
		if (save_result)
		{
			string line <- "" + total_population + "," + ave_walking_distance + "," + ave_walking_speed + "," + ave_travel_time + "," + egress_time;
			save (line) to: "../includes/paper/" + folder + system_result_filename type: "csv" header: false rewrite: false;
		}

		/* 
		loop ctr over: counters_list
		{
			if (length(zombie where (each.rcso_counter_no = ctr)) > 0)
			{
				write "Total Ped for Counter " + ctr + " :" + length(zombie where (each.rcso_counter_no = ctr));
				write "% Ped for Counter " + ctr + " :" + (length(zombie where (each.rcso_counter_no = ctr)) / total_zombie) * 100;
				write "RCSO AWD for Counter " + ctr + " :" + (sum((zombie where (each.rcso_counter_no = ctr)) collect (each.my_performance_record.walkingDistance)) / length(zombie where
				(each.rcso_counter_no = ctr))) with_precision 2;
				write "RCSO AWS for Counter " + ctr + " :" + (sum((zombie where (each.rcso_counter_no = ctr)) collect (each.my_performance_record.speed)) / length(zombie where
				(each.rcso_counter_no = ctr))) with_precision 2;
				write "RCSO ATT for Counter " + ctr + " :" + (sum((zombie where (each.rcso_counter_no = ctr)) collect (each.my_performance_record.travelTime)) / length(zombie where
				(each.rcso_counter_no = ctr))) with_precision 2;
			}

		}
		
*/

//	string line <- string(beta_cdf_alpha) + " " + beta_cdf_beta + " " + string(ave_walking_distance) + " " + ave_travel_time + " " + ave_walking_speed;
		//	save line to: "test.asc";

	}

	action update_colors
	{
		ask people
		{
			float alpha <- self.my_cell.density * 0.2;
			color <- rgb(255, 0, 0, alpha);
		}

	}

	action update_density
	{
		ask cell
		{
			density <- length(people inside self);
		}

	}

	//actions for performance computation
	action save_zombies
	{
	//write "saving zombies";
		ask zombie
		{
			string line <- string(time_expired) + "," + string(exit_cell);
			save (line) to: "../includes/paper/" + folder + "/results/" + length(zombie) + "_door_" + smethod + door + "_R" + run_no + "_evacuation_time.csv" rewrite: false type: "csv";
		}

	}

	/*action plot_egress{
		//clear the file
		save to: "../includes/paper/" + folder + "/" + length(zombie) + "_" + door + "_egress.csv" rewrite: true;
		
		float max_egress <- max(zombie collect (each.time_expired));
		write "max_egrees" + max_egress;
		int ctr <-0;
			loop c from: 0 to:  (max_egress*2) 
			{
			  float t <-  c * delta_time;
		      int zombie_total <- zombie count (each.time_expired = t);
		       ctr <- ctr + zombie_total;
		   
		       	string line <- string(t) + ","+ string(ctr);
		    	save (line) to: "../includes/paper/" + folder + "/" + length(zombie) + "_" + door + "_egress.csv" rewrite: false;
		
		    }
		
	}*/
	action plot_egress
	{
	//clear the file
	//	save to: "../includes/paper/" + folder + "/results/" + length(zombie)  + "_egress.csv" rewrite: false type:csv;
		float max_egress <- max(zombie collect (each.time_expired));
		loop list_ctr over: counters_list
		{
			int zombie_ctr <- 0;
			loop c from: 0 to: (max_egress + 1)
			{
				int zombie_total <- zombie count ((each.rcso_counter_no = list_ctr) and (each.time_expired <= c));
				string line <- string(c) + "," + string(zombie_total) + ",c" + list_ctr + "," + smethod + "," + door;
				save (line) to: "../includes/paper/" + folder + "/results/" + length(zombie) + "_egress.csv" rewrite: false type: "csv";
			}

		}

	}

} //end of global

//-------------start of species -------------------------------------
species people parent: people_base
{
	cell my_cell <- nil;
	cell my_source_cell <- nil;
	cell my_target_cell <- nil;
	cell my_previous_cell <- nil;
	bool is_dead <- false;
	int entry_angle <- 0;
	int exit_angle <- 0;
	list<list<float>> waypoints <- nil;
	list<list<float>> locations_in_cell <- nil;
	bool ready_to_move <- false;
	
	bool isInMesoZone <- false;
	bool can_move_out
	{
		if(isInMesoZone){
			return true;
		}
		float speed <- S(length(people inside my_cell) - 1);
		float next_position <- position_in_cell + (speed * delta_time);

		//meaning next position is outside the diameter
		if (next_position > cell_diameter)
		{
			return true;
		} else
		{
			return false;
		}

	}

	reflex update_life
	{
		if ((my_cell.basin in [1]))
		{
			is_dead <- true;
		}

	}

	list<cell> sort_by_p2enter (list<cell> cell_choices, matrix nav_matrix)
	{
		float current_spv <- float(nav_matrix[cell(self.location).grid_x, cell(self.location).grid_y]);
		loop g over: cell_choices
		{
			g.spv <- float(nav_matrix[g.grid_x, g.grid_y]);
		}

		list<float> spvs <- cell_choices collect (each.spv);
	
		loop g over: cell_choices
		{
		
			g.probability_to_enter <- ((I(g.density) * (1.0 - Navigation_bias)) + (N(g.spv, spvs) * Navigation_bias)) ;
		}

		return cell_choices sort_by (-each.probability_to_enter);
	}

	cell find_best_cell (list<cell> cell_choices, cell c)
	{

  
	//get nearest by distance to force agent to move lateral
		list<cell> top_cells <- cell_choices where (each.probability_to_enter = c.probability_to_enter);
		
		//return cell  closest to cell location of person
		//return top_cells with_min_of (each.cell_distance(cell(location)));
	  
	   cell top_spv <- top_cells with_max_of(each.spv);
	   
	   top_cells <- top_cells where (each.spv = top_spv.spv);
	     
	     
		return  top_cells with_min_of (each.cell_distance(cell(location)));
	}

	action mark_next_destination
	{
	//based on local
		if (!is_dead and can_move_out() and !just_in)
		{
			matrix nav_matrix <- nav_map[nav_id];
			list<cell> cell_choices <- nil;
			//cell choices include all neighbours
			cell_choices <- cell(location).neighbours;
			//sort in descending order
			cell_choices <- (sort_by_p2enter(cell_choices, nav_matrix));

			//my_target_cell <- find_best_destination(cell_choices, nav_matrix);
			//if target cell not full then mark as next destination
			if (stay_in_place)
			{
				my_target_cell <- find_best_cell(cell_choices, cell_choices[0]);
				if (my_target_cell.temp_cell_density < max_density)
				{
					my_target_cell.temp_cell_density <- my_target_cell.temp_cell_density + 1;
					ready_to_move <- true;
				} else
				{
					my_target_cell <- my_cell;
				}

			} else
			{
				loop c over: cell_choices
				{
					my_target_cell <- find_best_cell(cell_choices, c);
					if (my_target_cell.temp_cell_density < max_density)
					{
						my_target_cell.temp_cell_density <- my_target_cell.temp_cell_density + 1;
						ready_to_move <- true;
						break;
					}

				}

			}

			next_location <- my_target_cell.location;
		}
		//mark temp_cell
	}

	action move
	{

	//waypoint is time / x /y
		my_cell.density <- length(people inside (my_cell));
		float speed <- S(my_cell.density - 1);
		int angle <- 0;
		if (ready_to_move and !just_in)
		{

		//A,B,C A is midpoint B is vertical 12oclock c is moving angle
			angle <- abs(angle_between({ my_cell.grid_x, -my_cell.grid_y }, { my_cell.grid_x + 1, -my_cell.grid_y }, { my_target_cell.grid_x, -my_target_cell.grid_y }));

			//store angle moving out
			exit_angle <- angle;
			//add locations in cell to waypoints
			if (my_cell.is_trap)
			{
				loop g over: locations_in_cell
				{
					float r <- 0.0;
					if (g[1] <= cell_diameter / 2)
					{
						r <- (cell_diameter / 2 - g[1]) with_precision 2;
						float x <- getXCoordinate(my_cell.grid_x, entry_angle, r);
						float y <- getYCoordinate(-my_cell.grid_y, entry_angle, r);
						waypoints << [g[0], x, y];
					} else
					{
						r <- (g[1] - (cell_diameter / 2)) with_precision 2;
						float x <- getXCoordinate(my_cell.grid_x, exit_angle, r);
						float y <- getYCoordinate(-my_cell.grid_y, exit_angle, r);
						waypoints << [g[0], x, y];
					}

				}

			}
			//clear locations in cell
			locations_in_cell <- nil;
			my_previous_cell <- my_cell;
			my_cell <- my_target_cell;
			location <- my_cell.location;
			//cell_time_in <- (cycle+1) * delta_time;
			ready_to_move <- false;
			if (my_cell.is_counter)
			{
				rcso_counter_no <- my_cell.rcso_counter_no;
			}

			position_in_cell <- (position_in_cell + (speed * delta_time)) with_precision 2;
			position_in_cell <- (position_in_cell - cell_diameter) with_precision 2;
			if (position_in_cell > cell_diameter / 2)
			{
			//add midpoint	if no midpoint exist
				float new_d <- (position_in_cell - (cell_diameter / 2));
				if (locations_in_cell contains [(((cycle + 1) * delta_time) - (new_d / speed)) with_precision 2, cell_diameter / 2])
				{
				} else
				{
					locations_in_cell << [(((cycle + 1) * delta_time) - (new_d / speed)) with_precision 2, cell_diameter / 2];
				}

			}

			locations_in_cell << [(cycle + 1) * delta_time, position_in_cell];
			entry_angle <- abs(angle_between({ my_cell.grid_x, -my_cell.grid_y }, { my_cell.grid_x + 1, -my_cell.grid_y }, { my_previous_cell.grid_x, -my_previous_cell.grid_y }));
		} else if (!just_in)
		{
			position_in_cell <- (position_in_cell + (speed * delta_time)) with_precision 2;
			if (position_in_cell > cell_diameter / 2)
			{
			//check if midpoint exist if not add the midpoint in locations in cell
				float last_loc <- (last(locations_in_cell))[1];
				float last_time <- (last(locations_in_cell))[0];
				if (last_loc < cell_diameter / 2)
				{
					float new_d <- ((cell_diameter / 2) - last_loc) with_precision 2;
					locations_in_cell << [last_time + (new_d / speed) with_precision 1, cell_diameter / 2];
				}

			}

			locations_in_cell << [(cycle + 1) * delta_time, position_in_cell];
		}

		if (just_in)
		{
			just_in <- !just_in;
		}

	}



} //end of people


grid meso_cell file: meso_layer neighbors: 8 schedules: [] use_individual_shapes: false use_regular_agents: false use_neighbors_cache: false{
	rgb color;
	int meso_number;
	
	aspect base{
		
		if (grid_value > 0)
		{
		//		color <- #gray;
		}
		
		
	}
}

grid cell file: micro_layer neighbors: 8 schedules: [] use_individual_shapes: false use_regular_agents: false use_neighbors_cache: false
{
	int basin <- 0;
	int basin_group <- 0;
	int reward <- 0;
	int entry_angle <- 0;
	bool is_counter <- false;
	list<cell> neighbours <- self neighbors_at 1;
	list<geometry> walls;
	rgb color;
	int temp_cell_density <- 0;
	int density <- 0;
	float probability_to_enter <- 0.0;
	float spv <- 0.0;
	float average_population <- 0.0;
	int rcso_counter_no <- 0;

	//for performance
	bool is_trap <- false;
	bool is_swamp <- false;
	bool is_meso <- false;
	float cell_distance (cell c1)
	{
	
		float f <- topology(world) distance_between ([self.location, c1.location]);
		return f;
	}

	action update_color
	{
		if (is_trap)
		{
			//color <- #beige;
		   color <- #white;
		}

		if (basin = 1)
		{
			color <- # gray;
		}

		if (basin = 2)
		{
			color <- # lightgray;
		}

		if (is_counter)
		{
			color <- #white;
		}
		
		if (is_swamp){
			color <- #brown;
		}

	}

	aspect draw_walls
	{	
			loop g over: walls
			{
				
				draw (g + 0.1) color: # black;
			}		
	}
	//draw text: "" + grid_x + "," + grid_y color: # gray size: 0.3;
	//draw text: name color: # gray size: 0.3;
	
}

species meso_block{
	list<meso_cell> member_cells;
	point center;
	int density;
	
	
	aspect base{
		draw circle(20) at: center color:#red;
	}
	

}

grid spvcells file: micro_layer neighbors: 8 schedules: [] {
	float spv <- 0.0;
	list<geometry> walls;
	action update_color
	{
   
		
	}
	
	aspect draw_walls
	{	
		 //  if (cell[grid_x,grid_y].is_swamp){
		//	color <- #brown;
		//   }
		
			loop g over: walls
			{
				
				draw (g + 0.1) color: # black;
			}		
			draw  "" +cell[grid_x,grid_y].name color: # gray size: 0.3;
	}
	
}

experiment main type: gui until: cycle > 0
{
	parameter "Speed Denstiy Alpha" var: Sa min: 0.1 max: 10.0 step: 0.1 category: "Parameters : High alpha means high value is less distributed";
	parameter "Speed Denstiy Beta" var: Sb min: 0.1 max: 10.0 step: 0.1 category: "Parameters : High alpha means high value is less distributed";
	parameter "Interaction Alpha" var: Ia min: 0.1 max: 10.0 step: 0.1 category: "Parameters : High alpha means high value is less distributed";
	parameter "Interaction Beta" var: Ib min: 0.1 max: 10.0 step: 0.1 category: "Parameters : High alpha means high value is less distributed";
	parameter "Navigation Alpha" var: Na min: 0.1 max: 10.0 step: 0.1 category: "Parameters : High alpha means high value is less distributed";
	parameter "Navigation Beta" var: Nb min: 0.1 max: 10.0 step: 0.1 category: "Parameters : High alpha means high value is less distributed";
	parameter "Navigation Bias " var: Navigation_bias min: 0.1 max: 1.0 step: 0.1 category: "Parameters : High Navigation means high SPV preferred";
	parameter "flow" var: flow min: 1 max: 5;
	parameter "Folder" var: folder;
	parameter "Penalty" var: penalty;
	parameter "Population " var: total_population;
	parameter "Load Navigation" var: load_navigation;
	parameter "Save Zombies" var: save_zombies;
	parameter "Save NTXY " var: save_ntxy;
	parameter "Show result" var: showResult;
	parameter "search method" var: smethod;
	parameter "Stay In Place" var: stay_in_place;
	parameter "USe People Generator" var: use_people_generator;
	output
	{
		display main
		{
			grid cell;
			grid meso_cell;
			species meso_cell aspect:base;	
			species cell aspect: draw_walls;
		}
		
//		display spv{
//			grid spvcells lines:#black;
//			species spvcells aspect: draw_walls;
//		}

	}

}

experiment SA type: batch repeat: 5 keep_seed: false until: abs(actual_egress_time - ave_egress_time) = 0.0
{
	parameter "Speed Denstiy Alpha" var: Sa min: 0.1 max: 0.9 step: 0.1 category: "Parameters : High alpha means high value is less distributed";
	parameter "Speed Denstiy Beta" var: Sb min: 0.1 max: 0.9 step: 0.1 category: "Parameters : High alpha means high value is less distributed";
	parameter "Interaction Alpha" var: Ia min: 0.1 max: 0.9 step: 0.1 category: "Parameters : High alpha means high value is less distributed";
	parameter "Interaction Beta" var: Ib min: 0.1 max: 0.9 step: 0.1 category: "Parameters : High alpha means high value is less distributed";
	parameter "Navigation Alpha" var: Na min: 0.1 max: 0.9 step: 0.1 category: "Parameters : High alpha means high value is less distributed";
	parameter "Navigation Beta" var: Nb min: 0.1 max: 0.9 step: 0.1 category: "Parameters : High alpha means high value is less distributed";
	parameter "Navigation Bias " var: Navigation_bias min: 0.1 max: 0.9 step: 0.1 category: "Parameters : High Navigation means high SPV preferred";
	parameter "Egress to Search" var: actual_egress_time;
	parameter "Population" var: total_population;
	parameter "Load Navigation" var: load_navigation;
	method annealing temp_init: 1000 temp_end: 1 temp_decrease: 0.5 nb_iter_cst_temp: 50 minimize: abs(actual_egress_time - egress_time);
	float ave_egress_time <- 0.0;
	int ctr <- 5;
	float ec <- 0.0;
	action _step_
	{
		write
		"" + Ia + "," + Ib + "," + Sa + "," + Sb + "," + Na + "," + Nb + "," + Navigation_bias + "," + egress_time + "," + ave_travel_time + "," + ave_walking_distance + "," + ave_walking_speed;
		ec <- ec + egress_time;
		ctr <- ctr - 1;
		if (ctr = 0)
		{
		//get average
			ave_egress_time <- ec / 5;
			//write "average et: " + ave_egress_time;
			//reset
			ctr <- 5;
			ec <- 0.0;
		}

	}

}

experiment Tabu_Search type: batch repeat: 5 keep_seed: false until: abs(actual_egress_time - ave_egress_time) = 0.0
{
	parameter "Speed Denstiy Alpha" var: Sa min: 0.1 max: 0.9 step: 0.1 category: "Parameters : High alpha means high value is less distributed";
	parameter "Speed Denstiy Beta" var: Sb min: 0.1 max: 0.9 step: 0.1 category: "Parameters : High alpha means high value is less distributed";
	parameter "Interaction Alpha" var: Ia min: 0.1 max: 0.9 step: 0.1 category: "Parameters : High alpha means high value is less distributed";
	parameter "Interaction Beta" var: Ib min: 0.1 max: 0.9 step: 0.1 category: "Parameters : High alpha means high value is less distributed";
	parameter "Navigation Alpha" var: Na min: 0.1 max: 0.9 step: 0.1 category: "Parameters : High alpha means high value is less distributed";
	parameter "Navigation Beta" var: Nb min: 0.1 max: 0.9 step: 0.1 category: "Parameters : High alpha means high value is less distributed";
	parameter "Navigation Bias " var: Navigation_bias min: 0.1 max: 0.9 step: 0.1 category: "Parameters : High Navigation means high SPV preferred";
	parameter "Egress to Search" var: actual_egress_time;
	parameter "Population" var: total_population;
	parameter "Folder" var: folder;
	parameter "Save NTXY " var: save_ntxy;
	parameter "Save Zombies " var: save_zombies;
	parameter "Load Navigation" var: load_navigation;
	method tabu iter_max: 200 minimize: abs(actual_egress_time - egress_time);
	float ave_egress_time <- 0.0;
	int ctr <- 5;
	float ec <- 0.0;
	action _step_
	{
		write
		"" + Ia + "," + Ib + "," + Sa + "," + Sb + "," + Na + "," + Nb + "," + Navigation_bias + "," + egress_time + "," + ave_travel_time + "," + ave_walking_distance + "," + ave_walking_speed;
		ec <- ec + egress_time;
		ctr <- ctr - 1;
		if (ctr = 0)
		{
		//get average
			ave_egress_time <- ec / 5;
			//write "average et: " + ave_egress_time;
			//reset
			ctr <- 5;
			ec <- 0.0;
		}

	}

}

experiment Hill_climbing type: batch repeat: 5 keep_seed: false until: abs(actual_egress_time - ave_egress_time) = 0.0
{
	parameter "Speed Denstiy Alpha" var: Sa min: 0.1 max: 0.9 step: 0.1 category: "Parameters : High alpha means high value is less distributed";
	parameter "Speed Denstiy Beta" var: Sb min: 0.1 max: 0.9 step: 0.1 category: "Parameters : High alpha means high value is less distributed";
	parameter "Interaction Alpha" var: Ia min: 0.1 max: 0.9 step: 0.1 category: "Parameters : High alpha means high value is less distributed";
	parameter "Interaction Beta" var: Ib min: 0.1 max: 0.9 step: 0.1 category: "Parameters : High alpha means high value is less distributed";
	parameter "Navigation Alpha" var: Na min: 0.1 max: 0.9 step: 0.1 category: "Parameters : High alpha means high value is less distributed";
	parameter "Navigation Beta" var: Nb min: 0.1 max: 0.9 step: 0.1 category: "Parameters : High alpha means high value is less distributed";
	parameter "Navigation Bias " var: Navigation_bias min: 0.1 max: 0.9 step: 0.1 category: "Parameters : High Navigation means high SPV preferred";
	parameter "Egress to Search" var: actual_egress_time;
	parameter "Population" var: total_population;
	parameter "Folder" var: folder;
	parameter "Save NTXY " var: save_ntxy;
	parameter "Save Zombies " var: save_zombies;
	parameter "Load Navigation" var: load_navigation;
	method hill_climbing iter_max: 200 minimize: abs(actual_egress_time - egress_time);
	float ave_egress_time <- 0.0;
	int ctr <- 5;
	float ec <- 0.0;
	action _step_
	{
		write
		"" + Ia + "," + Ib + "," + Sa + "," + Sb + "," + Na + "," + Nb + "," + Navigation_bias + "," + egress_time + "," + ave_travel_time + "," + ave_walking_distance + "," + ave_walking_speed;
		ec <- ec + egress_time;
		ctr <- ctr - 1;
		if (ctr = 0)
		{
		//get average
			ave_egress_time <- ec / 5;
			//write "average et: " + ave_egress_time;
			//reset
			ctr <- 5;
			ec <- 0.0;
		}

	}

}

experiment Reactive_Tabu_Search type: batch repeat: 5 keep_seed: false until: abs(actual_egress_time - ave_egress_time) = 0.0
{
	parameter "Speed Denstiy Alpha" var: Sa min: 0.1 max: 0.9 step: 0.1 category: "Parameters : High alpha means high value is less distributed";
	parameter "Speed Denstiy Beta" var: Sb min: 0.1 max: 0.9 step: 0.1 category: "Parameters : High alpha means high value is less distributed";
	parameter "Interaction Alpha" var: Ia min: 0.1 max: 0.9 step: 0.1 category: "Parameters : High alpha means high value is less distributed";
	parameter "Interaction Beta" var: Ib min: 0.1 max: 0.9 step: 0.1 category: "Parameters : High alpha means high value is less distributed";
	parameter "Navigation Alpha" var: Na min: 0.1 max: 0.9 step: 0.1 category: "Parameters : High alpha means high value is less distributed";
	parameter "Navigation Beta" var: Nb min: 0.1 max: 0.9 step: 0.1 category: "Parameters : High alpha means high value is less distributed";
	parameter "Navigation Bias " var: Navigation_bias min: 0.1 max: 0.9 step: 0.1 category: "Parameters : High Navigation means high SPV preferred";
	parameter "Egress to Search" var: actual_egress_time;
	parameter "Population" var: total_population;
	parameter "Folder" var: folder;
	parameter "Save NTXY " var: save_ntxy;
	parameter "Save Zombies " var: save_zombies;
	parameter "Load Navigation" var: load_navigation;
	method reactive_tabu iter_max: 200 minimize: abs(actual_egress_time - egress_time);
	float ave_egress_time <- 0.0;
	int ctr <- 5;
	float ec <- 0.0;
	action _step_
	{
		write
		"" + Ia + "," + Ib + "," + Sa + "," + Sb + "," + Na + "," + Nb + "," + Navigation_bias + "," + egress_time + "," + ave_travel_time + "," + ave_walking_distance + "," + ave_walking_speed;
		ec <- ec + egress_time;
		ctr <- ctr - 1;
		if (ctr = 0)
		{
		//get average
			ave_egress_time <- ec / 5;
			//write "average et: " + ave_egress_time;
			//reset
			ctr <- 5;
			ec <- 0.0;
		}

	}

}

experiment GA type: batch keep_seed: true until: empty(people)
{
	parameter "Speed Denstiy Alpha" var: Sa min: 0.1 max: 0.9 step: 0.1 category: "Parameters : High alpha means high value is less distributed";
	parameter "Speed Denstiy Beta" var: Sb min: 0.1 max: 0.9 step: 0.1 category: "Parameters : High alpha means high value is less distributed";
	parameter "Interaction Alpha" var: Ia min: 0.1 max: 0.9 step: 0.1 category: "Parameters : High alpha means high value is less distributed";
	parameter "Interaction Beta" var: Ib min: 0.1 max: 0.9 step: 0.1 category: "Parameters : High alpha means high value is less distributed";
	parameter "Navigation Alpha" var: Na min: 0.1 max: 0.9 step: 0.1 category: "Parameters : High alpha means high value is less distributed";
	parameter "Navigation Beta" var: Nb min: 0.1 max: 0.9 step: 0.1 category: "Parameters : High alpha means high value is less distributed";
	//parameter "Denstiy Alpha" var: Da min: 0.1 max: 0.9 step: 0.1 category: "Parameters : High alpha means high value is less distributed";
	parameter "Navigation Bias " var: Navigation_bias min: 0.1 max: 0.9 step: 0.1 category: "Parameters : High Navigation means high SPV preferred";
	parameter "Population " var: total_population;
	parameter "Egress to Search" var: actual_egress_time;
	parameter "Load Navigation" var: load_navigation;
	method genetic minimize: ave_walking_speed max_gen: 200;
	float ave_egress_time <- 0.0;
	int ctr <- 5;
	float ec <- 0.0;
	init
	{
		showResult <- true;
	}

}

experiment FD type: batch repeat: 5 until: empty(people)
{
	parameter "Population" var: total_population min: 100 max: 10000 step: 100;
	//parameter "Save NTXY " var: save_ntxy;
	
	init
	{
		Navigation_bias <- 0.5;
		load_navigation <-true;
		
		
		
		
			write "rewriting csv";
			string line <- "pop,awd,aws,att,et";
			save (line) to: "../includes/paper/" + folder + system_result_filename type: "csv" header: false rewrite: true;
		

	}

}

