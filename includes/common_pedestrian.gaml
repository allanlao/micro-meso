/**
 *  commonpedestrian
 *  Author: allanlao
 *  Description: 
 */
model commonpedestrian


global
{
	
	
	float learning_rate <- 0.8;
	int max_iteration <- 200000;
	map<int, matrix> nav_map <- [];
	bool use_people_generator <- false;
	float Sa <- 1.0;
	float Sb <- 1.0;
	float Na <- 1.0;
	float Nb <- 1.0;
	float Navigation_bias <- 0.2;
	float Ia <- 1.0;
	float Ib <- 1.0;
	float cell_diameter <- 1.0;
	int max_density <- 5;
	float delta_time <- 0.5;
	float max_speed <- 1.4;

	//performance computation
	//performance indicators
	float ave_walking_distance <- 0.0;
	float ave_straight_distance_ratio <- 0.0;
	float ave_individual_delay <- 0.0;
	float ave_walking_speed <- 0.0;
	float ave_travel_time <- 0.0;
	list<performance> performance_list <- nil;
	list<float, int> scenario_result <- nil;
	int run_no <- 0;
	
	map basin_map_probability (int n, int k, float p)
	{
		list nlist <- nil;
		loop times: k
		{
			nlist << binomial(n, p);
		}

		return nlist frequency_of (each);
	}

}

species performance_record
{
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

species people_base
{
	rgb color <- # red;
	float time_expired <- 0.0;
	//float cell_time_in <- 0.0;
	float system_in <- 0.0;
	float system_out <- 0.0;
	int nav_id;
	point next_location <- nil;
	list<int> destination_list;
	int rcso_counter_no <- 0;
	float position_in_cell <- 0.0;
	bool just_in <- true;
	float getXCoordinate (int c, int angle, float radius)
	{
		float x <- c + radius * cos_rad(angle * # pi / 180);
		return x with_precision 2;
	}

	float getYCoordinate (int c, int angle, float radius)
	{
		float y <- c + radius * sin_rad(angle * # pi / 180);
		return y with_precision 2;
	}

	float time_start <- 0.0;
	float betaCDF (float x, float a, float b)
	{
		float value <- (incomplete_beta(a, b, x) / incomplete_beta(a, b, 1));
		return value with_precision 2;
	}

	float S (int density)
	{
		float x <- density / max_density;
		if (x > 1.0 or x < 0.0)
		{
			write "density is " + density + ", max_density: " + max_density;
		}

		float s <- max_speed * (1 - betaCDF(x, Sa, Sb)) with_precision 2;
		if (s <= 0.0)
		{
			s <- 0.1;
		}

		return s;
	}


	float N (float spv, list<float> spvs)
	{
		float min_spv <-  min(spvs collect(each));
		float max_spv <-  max(spvs collect(each));
		
		float z <- (spv - min_spv)/(max_spv-min_spv) ;
		return betaCDF(z, Na, Nb) with_precision 4;
	}

	float I(int density, list<int> densities)
	{
		int a <-  min(densities collect(each));
		int b <-  max(densities collect(each));
		
		float x <- (density - a)/(b-a) ;
		
		
	
		float value <- 1 - betaCDF(x, Ia, Ib);
		return value with_precision 2;
	}

}

species zombie
{
	float time_expired <- 0.0;
	int cycles_expired <- 0;
	list<list<float>> waypoints <- nil;
	int rcso_counter_no <- 0;
	performance_record my_performance_record <- nil;
	float terminal_time <- 0.0;
	string peoplename <- nil;
	string exit_cell <- nil;
}

species performance
{
	string agent_name <- nil;
	float time_start <- 0.0;
	float time_end <- 0.0;
}


