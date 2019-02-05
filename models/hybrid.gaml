model hybrid

global {
	string folder <- 'san-juan-meso';
	string system_result_filename <- 'system_result';
	file micro_file <- file("../includes/paper/" + folder + "/micro_layer.asc");
	file meso_file <- file("../includes/paper/" + folder + "/meso_layer.asc");
	file basin_file <- csv_file("../includes/paper/" + folder + "/basin.csv", ",");
	file trap_file <- file("../includes/paper/" + folder + "/trap.csv");
	file counter_file <- file("../includes/paper/" + folder + "/counters.csv");
	int population <- 100 min: 1 max: 1000;
	bool use_people_generator <- false;
	int people_ctr <- 0;
	int flow <- 3;
	string run_desc <- nil;
	bool load_navigation <- true;
	bool batch_mode <- false;
	list<cell> free_space <- nil;
	list<int> meso_group <- nil;
	list<cell> source_cell <- nil;
	list<int> counters_list <- nil;
	map prob_dist_map <- nil;
	int dist_k <- 200;
	int dist_n <- 0;
	float dist_p <- 0.5;
	int max_density <- 4;
	float base_circle_size <- 0.15;
	float size_inc <- 0.1;
	list<int> basin_group <- nil;
	float learning_rate <- 0.8;
	int max_iteration <- 200000;
	float max_speed <- 1.4;
	float cell_size <- 0.5;

	//float micro_cell_diameter <- (((cell_size * sqrt(2)) + cell_size)/2) with_precision 2;
	//float meso_cell_diameter <-  (((cell_size * 2 * sqrt(2)) + cell_size * 2)/2) with_precision 2;
	float micro_cell_diameter <- cell_size;
	float meso_cell_diameter <- cell_size * 2;
	//float cell_unit_size <- 0.5; //the actual size of each cell used for computation of distance.
	//time to travel 0.5 m/s	
	float delta_time <- (micro_cell_diameter / max_speed) with_precision 2;

	/*float Ia <- 0.8;
	float Ib <- 0.5; //higher value means lesser probility for denser cells
	float Sa <- 1.0;
	float Sb <- 0.1;
	float Na <- 0.7;
	float Nb <- 1.0;
	float Navigation_bias <- 0.7;
	*/
	float Ia <- 1.0;
	float Ib <- 1.0; //higher value means lesser probility for denser cells
	float Sa <- 1.0;
	float Sb <- 1.0;
	float Na <- 0.1;
	float Nb <- 1.0;
	float Navigation_bias <- 0.5;

	//performance indicators
	float ave_walking_distance <- 0.0;
	float ave_straight_distance_ratio <- 0.0;
	float ave_individual_delay <- 0.0;
	float ave_walking_speed <- 0.0;
	float ave_travel_time <- 0.0;
	float egress_time <- 0.0;
	bool meso_mode <- true;
	bool stay_in_place <- false;
	bool save_ntxy <- false;
	bool save_result <- true;
	bool save_performance <- true;
	bool save_rcso <- false;
	bool print_result <- true;
	geometry shape <- envelope(micro_file);

	init {
		write ("delata time is" + delta_time);

		//*******build the space ************************
		do build_basins;
		if (counter_file != nil) {
			do build_rcso_counters;
		}

		do build_walls;
		do build_trap;
		
		if (load_navigation) {
			do loadNavigation;
		} else {
			do build_navigation;
		}

		
		ask cell {
			if (meso_cell[grid_x, grid_y].grid_value > 0) {
				if (meso_mode) {
					is_meso_cell <- true;
				}

				group_no <- int(meso_cell[grid_x, grid_y].grid_value);
			}

		}

		meso_group <- remove_duplicates(list<int>(meso_cell collect (each.grid_value)));
		remove item: 0 from: meso_group;
		loop g over: meso_group {
			if (g > 0) {
				create meso_block {
					group_no <- g;
					member_cells <- cell where (each.group_no = g);
					center <- point((member_cells[0].location.x + member_cells[3].location.x) / 2, (member_cells[0].location.y + member_cells[3].location.y) / 2);
				}

			}

		}
      
		//determine the meso block where the cell belongs
		
		ask (cell where (each.group_no >0)) {
			block <- meso_block(one_of(meso_block where (each.group_no = group_no)));
		}
   
		free_space <- cell where (length(each.neighbours) > 1);
		ask cell {
			neighborhood_choices <- neighbours;
			if (is_meso_cell) {
				loop g over: block.member_cells {
					add all: g.neighbours to: neighborhood_choices;
				}

				neighborhood_choices <- remove_duplicates(neighborhood_choices);
				//remove member cells from choice; meso should be represented by 1 virtual cell
				neighborhood_choices <- (neighborhood_choices - block.member_cells) + self;
			}

		}

		//get angle
		//source cells are basin = type 2	 
		source_cell <- cell where (each.basin = 2);
		source_cell <- source_cell sort (each.basin);

		//add distribution to source cells;
		dist_n <- length(source_cell) - 1;
		prob_dist_map <- basin_map_probability(dist_n, dist_k, dist_p);
		
		if (!use_people_generator) {
		
			source_cell <- (cell - (cell where (each.basin = 1))) where (each.is_trap);
			create people number: population {
			//   write "creating" +name;
				my_navigation <- one_of(navigation);
				//this is the problem
				//my_source_cell <- one_of(source_cell where(each.getDensity() < each.getMaxDensity() ));
				bool f <- false;
				loop while: (!f) {
					cell tcell <- one_of(source_cell);
					if (tcell.getDensity() < tcell.getMaxDensity()) {
						my_source_cell <- tcell;
						f <- !f;
					}

				}

				if (my_source_cell.getDensity() >= my_source_cell.getMaxDensity()) {
					source_cell <- source_cell - my_source_cell;
				}
				// my_source_cell <- cell(30);
				my_cell <- my_source_cell;
				do moveTo(my_cell);

				//position in cell is the linear pt location
				//locations_in_cell list of position in cell
				if (my_cell.is_meso_cell) {
				//initial position is at center that is 0.5 of 1 meter
					position_in_cell <- (meso_cell_diameter / 2);
				} else {
					position_in_cell <- (micro_cell_diameter / 2);
				}
				//if people is created at once then initial position should change ;
				locations_in_cell << [(cycle) * delta_time, position_in_cell];
				entry_angle <- 0;
			} //create people
			
		}

	} //end of init


	//additional global variables
	map basin_map_probability (int n, int k, float p) {
		list nlist <- nil;
		loop times: k {
			nlist << binomial(n, p);
		}

		return nlist frequency_of (each);
	}

	//***************************************// 
	//my actions
	reflex create_people when: use_people_generator and (people_ctr < population) {
	//run generator
		loop times: flow {
			loop x over: prob_dist_map.keys {
				x <- int(x);
				if (length(people inside source_cell[x]) < source_cell[x].getMaxDensity()) {
					bool head <- flip((float(prob_dist_map[x]) / dist_k * 50 + 50) / 100);
					if (head and people_ctr < population) {
						create people number: 1 {
							my_source_cell <- source_cell[x];
							my_navigation <- one_of(navigation where (each.basin_group = 1)); //should be from destination list
							my_cell <- my_source_cell;
							do moveTo(my_cell);

							//position in cell is the linear pt location
							//locations_in_cell list of position in cell
							if (my_cell.is_meso_cell) {
							//initial position is at center that is 0.5 of 1 meter
								position_in_cell <- (meso_cell_diameter / 2);
							} else {
								position_in_cell <- (micro_cell_diameter / 2);
							}
							//if people is created at once then initial position should change ;
							locations_in_cell << [(cycle) * delta_time, position_in_cell];
							entry_angle <- 0;
						}

						people_ctr <- people_ctr + 1;
					}

				}

			}

		}

	}

	action build_rcso_counters {
	//cx1,y1,x2,y2,counter_no
		matrix m <- matrix(counter_file);
	
		if (m.rows <= 1) {
			return;
		}

		loop row_ctr from: 0 to: m.rows - 1 {
			int x1 <- 0;
			int x2 <- 0;
			int y1 <- 0;
			int y2 <- 0;
			x1 <- int(m[0, row_ctr]);
			y1 <- int(m[1, row_ctr]);
			x2 <- int(m[2, row_ctr]);
			y2 <- int(m[3, row_ctr]);
			loop col from: x1 to: x2 {
				loop row from: y1 to: y2 {
					cell[col, row].rcso_counter_no <- int(m[4, row_ctr]);
					cell[col, row].color <- #orange;
					cell[col, row].is_counter <- true;
				}

			}

			add int(m[4, row_ctr]) to: counters_list;
		}

	}

	action save_rcso_table {
		write "saving rcso to " + "../includes/paper/" + folder + "/results/rcso_" + population + "_" + time;
		int totalCtrs <- (cell where(each.is_counter));
		int rc1 <- 0;
		int rc2 <- 0;
		int ctr <- 0;
		ask zombie {
			if (rcso_counter_no = 1) {
				rc1 <- rc1 + 1;
			} else {
				rc2 <- rc2 + 1;
			}

			ctr <- ctr + 1;
			string line <- "" + (ctr) + "," + rc1 + "," + rc2;
			save (line) to: "../includes/paper/" + folder + "/results/rcso_" + population + "_" + run_desc type: "csv" header: false rewrite: false;
		}

	}

	action compute_performance {
		ask zombie {
			my_performance_record.travelTime <- 0.0;
			//ind walking distance
			loop ctr from: 0 to: length(waypoints) - 2 {
			//waypoint TXYP
			// chebyshev
				float d <- max([abs(waypoints[ctr][1] - waypoints[ctr + 1][1]), abs(waypoints[ctr][2] - waypoints[ctr + 1][2])]);

				//euclidean

				//	 float d <- sqrt((waypoints[ctr+1][1] - waypoints[ctr][1])^2 + (waypoints[ctr+1][2] - waypoints[ctr][2])^2);
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
		if (save_result) {
			string line <- "" + population + "," + ave_walking_distance + "," + ave_walking_speed + "," + ave_travel_time + "," + egress_time;
			save (line) to: "../includes/paper/" + folder + system_result_filename type: "csv" header: false rewrite: false;
		}

		loop ctr over: counters_list {
			if (length(zombie where (each.rcso_counter_no = ctr)) > 0) {
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

		//	string line <- string(beta_cdf_alpha) + " " + beta_cdf_beta + " " + string(ave_walking_distance) + " " + ave_travel_time + " " + ave_walking_speed;
		//	save line to: "test.asc";

	}

	action show_result {
		write "Population:" + population;
		write "parameters: NavBias>" + Navigation_bias;
		write "egress :" + egress_time + " att:" + ave_travel_time + " atd:" + ave_walking_distance + " ave walking speed:" + ave_walking_speed;
	}

	action build_trap {
		int x1 <- 0;
		int x2 <- 0;
		int y1 <- 0;
		int y2 <- 0;
		matrix m <- matrix(trap_file);
		x1 <- int(m[0, 0]);
		y1 <- int(m[1, 0]);
		x2 <- int(m[2, 0]);
		y2 <- int(m[3, 0]);
		loop c from: x1 to: x2 {
			loop r from: y1 to: y2 {
				cell[c, r].is_trap <- true;
			}

		}

	}

	action build_walls {
		write "building walls";
		ask cell {

		//top
			if (grid_value >= 8.0) {
				walls << line(shape.points where (each.y < location.y));
				neighbours <- neighbours where (each.location.y >= (location.y));
			}
			//left
			if (not even(grid_value)) {
				walls << line(shape.points where (each.x < location.x));
				neighbours <- neighbours where (each.location.x >= (location.x));
			}
			//bottom
			if (grid_value = 2.0 or grid_value = 3.0 or grid_value = 6.0 or grid_value = 7.0 or grid_value = 10.0 or grid_value = 11.0 or grid_value = 14.0 or grid_value = 15.0) {
				walls << line(shape.points where (each.y > location.y));
				neighbours <- neighbours where (each.location.y <= (location.y));
			}
			//right
			if ((grid_value >= 4.0 and grid_value <= 7.0) or (grid_value >= 12.0 and grid_value <= 15.0)) {
				walls << line(shape.points where (each.x > location.x));
				neighbours <- neighbours where (each.location.x <= (location.x));
			}

			add self to: neighbours at: 0;
		} //end

		ask cell {
			loop c over: copy(neighbours) {
				if not (self in c.neighbours) {
					using topology(world) {
						neighbours <- neighbours where ((each.location distance_to c.location) >= (each.location distance_to location));
					}

				}

			}

		}

		write "done";
	}

	action build_basins {
	//col,row,group,type,reward
	//type 1 = sink
	//type 2 = source
		matrix m <- matrix(basin_file);

		//skip header
		loop row_ctr from: 0 to: m.rows - 1 {
			int x1 <- 0;
			int x2 <- 0;
			int y1 <- 0;
			int y2 <- 0;
			x1 <- int(m[0, row_ctr]);
			y1 <- int(m[1, row_ctr]);
			x2 <- int(m[2, row_ctr]);
			y2 <- int(m[3, row_ctr]);
			loop col from: x1 to: x2 {
				loop row from: y1 to: y2 {
					cell[col, row].basin_group <- int(m[4, row_ctr]);
					cell[col, row].basin <- int(m[5, row_ctr]);
					cell[col, row].reward <- int(m[6, row_ctr]);
				}

			}

		}

	} //end of build basin
	action build_navigation {
	//sink cells = 1
	//basin_cells <- cell where (each.basin = 1);
		basin_group <- remove_duplicates(cell collect each.basin_group);
		basin_group <- basin_group where (each > 0);
		loop bg over: basin_group {
			create navigation {
				basin_group <- bg;
				nav_matrix <- generate_nav();
			}

		}

	}

	action loadNavigation {
		basin_group <- remove_duplicates(cell collect each.basin_group);
		basin_group <- basin_group where (each > 0);
		loop bg over: basin_group {
			file nav_file <- csv_file("../includes/paper/" + folder + "/navigation/" + bg + "_nav.csv", ",");
			matrix qm <- matrix(nav_file);
			create navigation {
				basin_group <- bg;
				nav_matrix <- qm;
			}

		}

	}

	action save_ntxy {
		string filename <- "/results/" + population + "_ntxyfile.csv";
		save ["n", "t", "x", "y", "p"] to: "../includes/paper/" + folder + filename type: "csv" header: false rewrite: true;
		ask zombie {
			loop ctr from: 0 to: length(waypoints) - 1 {
				string line <- name + "," + waypoints[ctr][0] + "," + waypoints[ctr][1] + "," + waypoints[ctr][2] + "," + population;
				save (line) to: "../includes/paper/" + folder + filename type: "csv" header: false rewrite: false;
			}

		}

	}

	action resolve_conflicts (list<cell> overflowing_micro, list<meso_block> overflowing_meso_blocks) {

	//check all micro cell first
		loop g over: overflowing_micro {
		//for every people not ready to move they have to fight it out with the place
			int total_safe <- g.tempPeople count (each.locked_position);
			list<people> fighting_people <- g.tempPeople where (!each.locked_position);

			//remove priority people from fighting people

			//fighting_people <- fighting_people - (fighting_people where(each.my_cell = each.my_target_cell)); 
			int slot <- g.getMaxDensity() - total_safe;

			//pick <<slot>> from among the fighting people
			list<people> winning_people <- nil;
			if (slot > 0) {
				winning_people <- slot among fighting_people;
			}

			//losing people to set my_target_to_nil and decrement tempdenstiy for the cell
			ask (fighting_people - winning_people) {
				ask my_target_cell {
					do remTempPeople(myself);
				}

				if (stay_in_place) {
					my_target_cell <- my_cell;
				} else {
					my_cell_choices <- my_cell_choices - my_target_cell;

					//move to next best location
					my_target_cell <- my_cell_choices[0];
				}

				if (my_target_cell = my_cell) {
					locked_position <- true;
				}

				ask my_target_cell {
					do addTempPeople(myself);
				}

			}

		} //end of loop g

		//then all the meso blocks;
		loop g over: overflowing_meso_blocks {
			list<people> fighting_people <- g.tempPeople where (!each.locked_position);
			int total_safe <- g.tempPeople count (each.locked_position);

			//for every people not ready to move they have to fight it out with the place
			int slot <- g.getMaxDensity() - total_safe;
			//pick <<slot>> from among the fighting people
			list<people> winning_people <- nil;
			if (slot > 0) {
				winning_people <- slot among fighting_people;
			}

			ask (fighting_people - winning_people) {
				ask my_target_cell {
					do remTempPeople(myself);
				}

				if (stay_in_place) {
					my_target_cell <- my_cell;
				} else {
					my_cell_choices <- my_cell_choices - my_target_cell;

					//move to next best location
					my_target_cell <- my_cell_choices[0];
				}

				if (my_target_cell = my_cell) {
					locked_position <- true;
				}

				ask my_target_cell {
					do addTempPeople(myself);
				}

			}

		} //overflowing meso;
	}

	//******** end of actions ***************//

	//**global reflexes ****************************/ 
	reflex set_density {
		ask cell {
			current_density <- getDensity();
		}

	}

	reflex remove_dead_people_in_micro_cell {
		ask people where (each.is_dead) {
			//if (my_cell.is_trap) {
				waypoints << [cycle * delta_time, my_cell.location.x, -my_cell.location.y];
			//}

			create zombie {
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

	reflex remove_dead_people_in_meso_cell {
		ask people where (each.my_cell.is_meso_cell) {
			if (is_dead and canMoveOut) {
				float speed <- S(my_cell.current_density - 1, my_cell.getMaxDensity());
				//		position_in_cell <- (cell_diameter) with_precision 2;
				if (position_in_cell > meso_cell_diameter / 2) {
				//check if midpoint exist if not add the midpoint in locations in cell
					float last_loc <- (last(locations_in_cell))[1];
					float last_time <- (last(locations_in_cell))[0];
					if (last_loc < meso_cell_diameter / 2) {
						float new_d <- ((meso_cell_diameter / 2) - last_loc) with_precision 2;
						locations_in_cell << [last_time + ((new_d / speed) with_precision 2), meso_cell_diameter / 2];
					}

				}
				//not including is trap because i want all areas to be counted
			//	if (my_cell.is_trap and my_cell.is_meso_cell) {
				if (my_cell.is_meso_cell) {
					locations_in_cell << [(cycle + 1) * delta_time, position_in_cell];
					loop g over: locations_in_cell {
						float r <- 0.0;
						//g[1] is the x axis
						//is x axis pt is greater that mid of diameter then its going out
						//x,y is a point in the circle inside the cell
						if (g[1] <= meso_cell_diameter / 2) {
							r <- (meso_cell_diameter / 2 - g[1]) with_precision 2;
							float x <- getXCoordinate(my_cell.getCenter().x, entry_angle, r);
							float y <- getYCoordinate(-my_cell.getCenter().y, entry_angle, r);
							waypoints << [g[0], x, y];
						} else {
							r <- (g[1] - meso_cell_diameter / 2) with_precision 2;
							float x <- getXCoordinate(my_cell.getCenter().x, exit_angle, r);
							float y <- getYCoordinate(-my_cell.getCenter().y, exit_angle, r);
							waypoints << [g[0], x, y];
						}

					}

				}

				create zombie {
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

	//REFLEX NO 1. FOR MOVEMENT
	// current density prior to any movement. this is the density at the time a people is deciding where to go.
	reflex reset_people_and_cell_status {
		ask people {
			locked_position <- false;
			my_target_cell <- nil;
			canMoveOut <- can_move();
		}

		ask cell {
			do setTempPeople(nil);
			current_density <- getDensity();
		}

	}

	reflex update_people_who_cannot_move {
		
		ask people where (!each.canMoveOut) {
		//just make them ready to move
			my_target_cell <- my_cell;
			//there are not moving so locked them out
			locked_position <- true;
			ask my_target_cell {
				do addTempPeople(myself);
			}

		}

	}

	reflex update_people_who_can_move {
	
		list<people> people_who_can_move <- people where (each.canMoveOut);
		//first time people who can move are moved. this is done only once.
		
		ask people_who_can_move {
			my_cell_choices <- my_cell.neighborhood_choices;
			my_cell_choices <- (sort_by_p2enter(my_cell_choices, my_navigation.nav_matrix));
			my_target_cell <- my_cell_choices[0];
			if (my_target_cell = my_cell) {
				locked_position <- true;
			}

			ask my_target_cell {
				do addTempPeople(myself);
			}

		}

		int lctr <- 0;
		
		loop while: true {

		//do the conflict resolution here
			list<cell> overflowing_micro <- cell where (each.getTempPeopleDensity() > each.getMaxDensity() and !each.is_meso_cell);
			list<meso_block> overflowing_meso_blocks <- meso_block where (each.getTempPeopleDensity() > each.getMaxDensity());
			if ((length(overflowing_micro) + length(overflowing_meso_blocks)) > 0) {
				do resolve_conflicts(overflowing_micro, overflowing_meso_blocks);
			} else {
				break;
			}

		}
    
	}

	reflex move_to_next_destination {
		
		list<people> people_can_move <- people where (each.canMoveOut);
		list<people> people_cannot_move <- people where (!each.canMoveOut);
		ask people_cannot_move {
			do update_waypoints(my_cell, my_target_cell);
			do move_in_place;
		}

		ask people_can_move {
			do update_waypoints(my_cell, my_target_cell);
			do move_agents;
		}
		

	}

	reflex terminate when: empty(people) {
		if (save_ntxy) {
		//	do save_ntxy;
		}

		if (save_performance) {
			do compute_performance;
		}

		if (save_rcso) {
		//	do save_rcso_table;
		}

		if (print_result) {
			do show_result;
		}

		if (!batch_mode) {
			do halt;
		}

	}
	//**end of global reflexes ****************************/ 

	//***************start of functions ***************//


	//*************end of functions *******************//

	//***********start of global actions*****************//

	//***********end of global actions************/
} //end of global
species navigation {
	file basin_file <- csv_file("../includes/paper/" + folder + "/basin.csv", ",");
	int basin_group;
	matrix nav_matrix;
	int penalty <- 0;
	int default_reward <- 0;
	float learning_rate <- 0.8;
	int max_iteration <- 200000;
	matrix generate_nav {
		int iter <- 0;
		int cols <- cell max_of (each.grid_x) + 1;
		int rows <- cell max_of (each.grid_y) + 1;
		list<list<float>> q_list <- [];
		//	list<cell> von_neighbours <- [];
		list<list<cell>> state_action_list <- [];
		matrix reward_matrix <- default_reward as_matrix ({cols, rows});
		matrix q_matrix <- 0.0 as_matrix ({cols, rows});

		//cells of the same group
		//	list<cell> basins <- cell where ((int(each.basin / 10)) = int(od_matrix[1, od_ctr])); //od_matrix col1(destination) row0
		list<cell> basins <- cell where (each.basin_group = basin_group);
		float current_sum <- 0.0;
		float prev_sum <- 0.0;
		int counter <- 0;

		//build reward

		//assign penalty to walls
		loop g over: cell {
			int p <- 0;
			if (length(g.neighbours) < 9) {
				p <- penalty;
			}

			if (g.is_swamp) {
				p <- p + penalty;
			}

			reward_matrix[g.grid_x, g.grid_y] <- -p;
		}

		//get reward of basins from rewards file
		loop g over: basins {
			reward_matrix[g.grid_x, g.grid_y] <- g.reward;
		}

		//build qmatrix
		//von neuman
		ask cell {

		//	add neighbours collect (reward_matrix[each.grid_x, each.grid_y]) to: q_list1;
			add neighbours collect (0.0) to: q_list;
			add neighbours to: state_action_list;
		}

		//iterate up to max_iteration  or until diffence between previous and next state is nil
		loop i from: 0 to: max_iteration {

		//state_action_list is the neighborhood
			loop r from: 0 to: length(state_action_list) - 1 {
				cell current_state <- state_action_list[r][0];
				//loop for von neuman
				loop c from: 0 to: length(state_action_list[r]) - 1 {
				//q(s,a) per action and determine the maximum
					cell next_state <- state_action_list[r][c];
					float cd <- current_state.location distance_to next_state.location;
					float reward <- float(reward_matrix[next_state.grid_x, next_state.grid_y]);
					if (cd > 1.0) {
						reward <- 0;
					}

					float max_q <- max(q_list[int(next_state)] collect (each)); //what is the max reward
					q_list[r][c] <- floor(reward + learning_rate * max_q);
				}

				//then loop for moore - von neuman

			} //end of i state action cycle
			loop g over: q_list {
				current_sum <- sum(g collect (each));
			}

			if (abs(prev_sum - current_sum) <= 0.000001) {
				counter <- counter + 1;
			} else {
				counter <- 0;
			}

			if (counter > 100) {
				write "breaking qlearning iteration";
				break;
			}

			prev_sum <- current_sum;
			iter <- i;
		} //max_iter


		//convert q_list to navigation matrix
		int ctr <- 0;

		//center
		loop row from: 0 to: q_matrix.rows - 1 {
			loop col from: 0 to: q_matrix.columns - 1 {
			//get center
				q_matrix[col, row] <- q_list[ctr][0];
				ctr <- ctr + 1;
			}

		}

		//max_spv <- max(q_matrix);

		//center
		loop row from: 0 to: q_matrix.rows - 1 {
			loop col from: 0 to: q_matrix.columns - 1 {
			//get center
				q_matrix[col, row] <- q_matrix[col, row];
				cell[col, row].spv_value <- q_matrix[col, row];
			}

		}

		do save_matrix(q_matrix, "../includes/paper/" + folder + "/navigation/" + basin_group + "_nav.csv");
		return q_matrix;
	}

	action save_matrix (matrix m, string f) {
		save m to: f rewrite: true type: "csv" header: "false";
	}

}

species people {
	rgb mycolor <- rgb(rnd(255), rnd(255), rnd(255));
	cell my_cell <- nil;
	cell my_source_cell <- nil;
	cell my_target_cell <- nil;
	cell my_previous_cell <- nil;
	bool is_dead <- false;
	int entry_angle <- 0;
	int exit_angle <- 0;
	list<list<float>> waypoints <- nil;
	list<list<float>> locations_in_cell <- nil;
	list<cell> my_cell_choices <- nil;
	bool locked_position <- false;
	bool canMoveOut <- false;
	navigation my_navigation;

	//point my_next_location <- nil;
	int rcso_counter_no <- 0;
	float position_in_cell <- 0.0;
	float random_number <- 0.0;

	//cell choice /neigbourhood will vary depending on state
	//micro is the same meso is union of neighbours of member cells
	bool can_move {
		if (!my_cell.is_meso_cell) {
			return true;
		} else {
			float speed <- S(my_cell.getDensity() - 1, my_cell.getMaxDensity());
			float next_position <- position_in_cell + (speed * delta_time);

			//meaning next position is outside the diameter
			if (next_position > meso_cell_diameter) {
				return true;
			} else {
				return false;
			}

		}

	}

	reflex update_life {
		if ((my_cell.basin in [1])) {
			is_dead <- true;
		}

	}

	action setLocation (point loc) {
		location <- my_cell.location;
	}

	float getXCoordinate (float c, int angle, float radius) {
		float x <- c + (radius * cos_rad(angle * #pi / 180));
		return x with_precision 2;
	}

	float getYCoordinate (float c, int angle, float radius) {
		float y <- c + (radius * sin_rad(angle * #pi / 180));
		return y with_precision 2;
	}

	float time_start <- 0.0;
	float betaCDF (float x, float a, float b) {
		float value <- (incomplete_beta(a, b, x) / incomplete_beta(a, b, 1));
		return value with_precision 2;
	}

	float S (int density, int md) {
		float x <- density / md;
		float s <- max_speed * (1 - betaCDF(x, Sa, Sb)) with_precision 2;
		if (s <= 0.0) {
			s <- 0.1;
		}

		return s;
	}

	float N (float spv, list<float> spvs) {
		float min_spv <- min(spvs collect (each));
		float max_spv <- max(spvs collect (each));
		float z;
		if (max_spv - min_spv = 0) {
			z <- 0.0;
		} else {
			z <- (spv - min_spv) / (max_spv - min_spv);
		}

		return betaCDF(z, Na, Nb) with_precision 2;
	}

	float I (int density, list densities) {
		int a <- min(densities collect (each));
		int b <- max(densities collect (each));
		float x;
		if (a = b) {
			x <- 0.0;
		} else {
			x <- (density - a) / (b - a);
		}

		float value <- 1 - betaCDF(x, Ia, Ib);
		return value with_precision 2;
	}

	//navigation matrix can be different among people depending on OD.
	//navigation is contained in nav.nav_matrix
	list<cell> sort_by_p2enter (list<cell> cell_choices, matrix nav_matrix) {
		float current_spv <- float(nav_matrix[cell(self.location).grid_x, cell(self.location).grid_y]);
		loop g over: cell_choices {
			g.spv_value <- float(nav_matrix[g.grid_x, g.grid_y]);
		}

		list<float> spvs <- cell_choices collect (each.getSpv());
		list<int> densities <- cell_choices collect (each.getDensity());
		loop g over: cell_choices {
			g.probability_to_enter <- (I(g.getDensity(), densities) * (1.0 - Navigation_bias)) + (N(g.getSpv(), spvs) * Navigation_bias);
		}

		return cell_choices sort_by (-each.probability_to_enter);
	}

	action moveTo (cell target) {
	//	if(target.getDensity() < target.getMaxDensity()){
		my_previous_cell <- my_cell;
		my_cell <- target;
		do setLocation(my_cell.location);

		//	}
	}

	action update_waypoints (cell current, cell target) {
		float diameter;
		if (current.is_meso_cell) {
			diameter <- meso_cell_diameter;
		} else {
			diameter <- micro_cell_diameter;
		}
		// if target same as current just add waypoint same place different time
		//position is the same

		//july 1,2-18 this has to be modified to stay in previous place not in the center
		if (current = target) {
			if (length(waypoints) > 0) {
				list<float> w <- last(waypoints);
				//put it in locations instead of waypoint..uncomment waypoint in this code wont work

				//	waypoints << [cycle*delta_time, w[1], w[2]];

			} else {
			//this will be the initial location of the agent inside a cell..its the center
			//		waypoints << [cycle*delta_time, current.getCenter().x, -current.getCenter().y];
			}

			return;
		}

		//waypoint is time / x /y
		float speed <- S(current.current_density - 1, current.getMaxDensity());
		int angle <- 0;

		//A,B,C A is midpoint B is vertical 12oclock c is moving angle
		exit_angle <-
		abs(angle_between({current.getCenter().x, -current.getCenter().y}, {current.getCenter().x + 1, -current.getCenter().y}, {target.getCenter().x, -target.getCenter().y}));
		//store angle moving out
		//people entering this stage contains locations_in_cell from previous cell (still the my cell not yet changed to targt cell)
		//add all locations in cell to waypoints
		//if (current.is_trap) {
			loop g over: locations_in_cell {

			//g[0] is time
			//g[1] is position in cell at that time
				float r <- 0.0;
				r <- abs(diameter / 2 - g[1]) with_precision 2;
				//if g[1] is <= halfway
				//radius is now full radius - distance travelled inside meso cell
				if (g[1] <= diameter / 2) {
					angle <- entry_angle;
				} else {
					angle <- exit_angle;
				}

				float x <- getXCoordinate(current.getCenter().x, angle, r);
				float y <- getYCoordinate(-current.getCenter().y, angle, r);
				waypoints << [g[0], x, y];
			} //after all locations in cell

			//determine the init poistion of people in the next cell
			//clear locations in cell
			locations_in_cell <- nil;
			position_in_cell <- (position_in_cell + (speed * delta_time)) with_precision 2;
			position_in_cell <- (position_in_cell - diameter) with_precision 2;
			locations_in_cell << [(cycle + 1) * delta_time, position_in_cell];
			entry_angle <-
			abs(angle_between({target.getCenter().x, -target.getCenter().y}, {target.getCenter().x + 1, -target.getCenter().y}, {current.getCenter().x, -current.getCenter().y}));
	//	} //end if trap

	} //end update waypoints
	action move_agents {
	//if current location of cell is meso
		do moveTo(my_target_cell);
		//update all
		//cell_time_in <- (cycle+1) * delta_time;
		if (my_cell.is_counter) {
			rcso_counter_no <- my_cell.rcso_counter_no;
		}

	} //end of move_agents
	action move_in_place {
		float speed <- S(my_cell.current_density - 1, my_cell.getMaxDensity());
		position_in_cell <- (position_in_cell + (speed * delta_time)) with_precision 2;

		//agent is past the midpoint
		//so check if there is an entry for midpoint already
		//make sure time isnt too close at precision 2 decimal
		if (position_in_cell > meso_cell_diameter / 2) {
		//check if midpoint exist if not add the midpoint in locations in cell
			float last_loc <- (last(locations_in_cell))[1];
			float last_time <- (last(locations_in_cell))[0];
			if (last_loc < meso_cell_diameter / 2) {
				float new_d <- ((meso_cell_diameter / 2) - last_loc) with_precision 2;
				float new_t <- (last_time + (new_d / speed)) with_precision 2;
				float tt <- ((cycle + 1) * delta_time) with_precision 2;
				if (tt != new_t) {
					locations_in_cell << [last_time + (new_d / speed) with_precision 2, meso_cell_diameter / 2];
				}

			}

		}

		locations_in_cell << [(cycle + 1) * delta_time, position_in_cell];
	}

	aspect default {
		if (!my_cell.is_meso_cell) {
			draw circle(base_circle_size + 0.05) color: #yellow;
			draw circle(base_circle_size) color: mycolor;
		} else {
			int x <- my_cell.getDensity();
			if (x > 0) {
				draw circle(base_circle_size + ((x - 1) * size_inc)) at: my_cell.getCenter() color: mycolor;
			}

		}

	}

} //end of people
species zombie {
	float time_expired <- 0.0;
	int cycles_expired <- 0;
	list<list<float>> waypoints <- nil;
	int rcso_counter_no <- 0;
	performance_record my_performance_record <- nil;
	float terminal_time <- 0.0;
	string peoplename <- nil;
	string exit_cell <- nil;
}

species performance {
	string agent_name <- nil;
	float time_start <- 0.0;
	float time_end <- 0.0;
}

species performance_record {
	list<float> insWalkingDistance <- nil;
	list<float> insWalkingDisplacement <- nil;
	list<float> insVelocity <- nil;
	list<float> insSpeed <- nil;
	float walkingDistance <- 0.0;
	float walkingDisplacement <- 0.0;
	float straightLineDistance <- 0.0;
	float straightDistanceRatio <- 0.0;
	float velocity <- 0.0;
	float speed <- 0.0;
	float walkingDelay <- 0.0;
	float walkingAngle <- 0.0;
	float movingDirection <- 0.0;
	float straightLineSpeed <- 0.0;
	float travelTime <- 0.0;
	float dt <- 0.0;
}

species meso_block {
	int group_no;
	list<cell> member_cells;
	point center;
	//int density <- 0 ;
	list<people> tempPeople <- nil;
	int getTempPeopleDensity {
	//return temp_cell_density;
		return length(tempPeople);
	}

	int getMaxDensity {
		return max_density;
	}

	aspect base {

	//location is x,y,z
		int x <- one_of(member_cells).getDensity();
		if (x > 0) {

		//draw circle(base_circle_size + ( (x-1) * size_inc)) at: center color:#gray;

		}

	}

}

grid meso_cell file: meso_file neighbors: 8 schedules: [] use_individual_shapes: false use_regular_agents: false use_neighbors_cache: false {
	rgb color;
	int density <- 0;

	aspect default {
		if (grid_value > 0) {
			color <- rgb(242, 242, 242);
		}

	}

}

grid cell file: micro_file neighbors: 8 schedules: [] use_individual_shapes: false use_regular_agents: false use_neighbors_cache: false {
	list<cell> neighbours <- self neighbors_at 1;
	list<geometry> walls;
	int group_no <- 0;
	meso_block block <- nil;
	float spv_value <- 0.0;
	int temp_cell_density <- 0;
	list<cell> neighborhood_choices;
	float probability_to_enter <- 0.0;
	int rcso_counter_no <- 0;
	int basin <- 0;
	int basin_group <- 0;
	int reward <- 0;
	int entry_angle <- 0;
	int current_density <- 0;
	bool is_swamp <- false;
	bool is_trap <- false;
	bool is_counter <- false;
	bool is_meso_cell <- false;
	bool is_free_space <- true;
	list<people> tempPeople <- nil;
	point getCenter {
		if (is_meso_cell) {
			return block.center;
		} else {
			return location;
		}

	}

	string getState {
		if (is_meso_cell) {
			return "meso_mode";
		} else {
			return "micro_mode";
		}

	}

	int getMaxDensity {
		if (is_meso_cell) {
			return max_density;
		} else {
			return 1;
		}

	}

	action setSpv (float v) {
		spv_value <- v;
	}

	float getSpv {
		if (is_meso_cell) {
			return block.member_cells max_of (each.spv_value);
		} else {
			return spv_value;
		}

	}

	action setTempPeople (list<people> p) {
		if (is_meso_cell) {
			block.tempPeople <- p;
		} else {
			tempPeople <- p;
		}

	}

	int getTempPeopleDensity {
		if (is_meso_cell) {
			return length(block.tempPeople);
		} else {
			return length(tempPeople);
		}

	}

	int getDensity {
		int x <- 0;
		if (is_meso_cell) {
			loop g over: block.member_cells {
				x <- x + length(people inside (g));
			}

		} else {
			x <- length(people inside (self));
		}

		return x;
	}

	/*action incDensity{
		if(is_meso_cell){
			block.density <- block.density + 1;
		}
	}*/
	action remTempPeople (people p) {
		if (is_meso_cell) {
			remove p from: block.tempPeople;
		} else {
			remove p from: tempPeople;
		}

	}

	action addTempPeople (people p) {
		if (is_meso_cell) {
			add p to: block.tempPeople;
		} else {
			add p to: tempPeople;
		}

	}

	float getMinSpv {
		return min(neighbours collect (each.getSpv()));
	}

	float cell_distance (cell c1) {
		float f <- topology(world) distance_between ([self.location, c1.location]);
		return f;
	}

	aspect default {
		if (group_no = 10) {
		// 	color <- #blue;
		}

		loop g over: walls {
			draw (g + 0.01) color: #black;
		}

		//	if (!is_meso_cell and getDensity() >0){
		//			draw circle(base_circle_size )  color:#gray;
		//		}
		if (basin = 1) {
			color <- #green;
		} else if (basin = 2) {
			color <- #blue;
		}

		if (is_trap) {
		//	color <- #beige;
		}

	}

}

experiment main type: gui {
	parameter "Stay in Place " var: stay_in_place;
	parameter "bathc mode " var: batch_mode <- false;
	parameter "run number " var: run_desc;
	parameter "Load Navigation" var: load_navigation <- false;
	//	parameter "Population " var: population;
	parameter "Population" var: population <- 100;
	output {
		display default {
			grid cell lines: #beige;
			//grid meso_cell lines:#beige;
			species meso_block aspect: base;
			species cell aspect: default;
			species people aspect: default;
			//species meso_cell;

		}

	}

}

experiment Tabu_Search type: batch repeat: 5 keep_seed: true until: empty(people) {
	parameter "Speed Denstiy Alpha" var: Sa min: 0.1 max: 1.0 step: 0.1 category: "Parameters : High alpha means high value is less distributed";
	parameter "Speed Denstiy Beta" var: Sb min: 0.1 max: 1.0 step: 0.1 category: "Parameters : High alpha means high value is less distributed";
	parameter "Interaction Alpha" var: Ia min: 0.1 max: 1.0 step: 0.1 category: "Parameters : High alpha means high value is less distributed";
	parameter "Interaction Beta" var: Ib min: 0.1 max: 1.0 step: 0.1 category: "Parameters : High alpha means high value is less distributed";
	parameter "Navigation Alpha" var: Na min: 0.1 max: 1.0 step: 0.1 category: "Parameters : High alpha means high value is less distributed";
	parameter "Navigation Beta" var: Nb min: 0.1 max: 1.0 step: 0.1 category: "Parameters : High alpha means high value is less distributed";
	parameter "Navigation Bias " var: Navigation_bias min: 0.3 max: 0.7 step: 0.1 category: "Parameters : High Navigation means high SPV preferred";
	parameter "Population" var: population <- 100;
	parameter "Save NTXY " var: save_ntxy <- false;
	parameter "Save RCSO " var: save_rcso <- false;
	parameter "Load Navigation" var: load_navigation <- true;
	parameter "batch Mode " var: batch_mode <- true;
	parameter "compute performance" var: save_performance <- true;
	method tabu iter_max: 500 tabu_list_size: 5 minimize: ave_walking_distance;

	reflex end_of_runs {
		string
		line <- "" + Ia + "," + Ib + "," + Sa + "," + Sb + "," + Na + "," + Nb + "," + Navigation_bias + "," + egress_time + "," + ave_travel_time + "," + ave_walking_distance + "," + ave_walking_speed;
		save (line) to: "../includes/paper/" + folder + "/parameters/atd_" + population type: "csv" header: false rewrite: false;
		write line;
	}

}

experiment repeat_5_times type: batch repeat: 5 keep_seed: true until: empty(people) {

	reflex end_of_runs {
		ask simulation {
			write "sample zombie";
			write sample(zombie);
		}

	}

} 

