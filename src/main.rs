use rand::Rng;
use std::fs::OpenOptions;
use std::io::{BufWriter, Write};
use std::fs::File;
use rand_distr::Distribution;
const BINDINGCHANCE: f64 = 0.5;
//const RSQRD: f64 = 0.05;
const RSQRD: f64 = 0.0025; // 0.05^2

fn main() -> Result<(), Box<dyn std::error::Error>>{
    //Ready the writing part of the code.
    let file = OpenOptions::new()
        .append(true)  // Enable appending
        .create(true)  // Create the file if it doesn't exist
        .open("BindingTarget.txt")?;

    let mut writer = BufWriter::new(file);  

    let file2 = OpenOptions::new()
    .append(true)  // Enable appending
    .create(true)  // Create the file if it doesn't exist
    .open("TTD.txt")?;
    let mut writer2 = BufWriter::new(file2);  


    let mut boundpercent:Vec<i32> = vec![0;50000];
    let number:usize = 1000;
    let mut boundto = vec![-1;number];
    let mut done:usize = 0;
    let dna = placedna(1);
    let mut particles = Placeparticles(number, &mut dna.clone());
    let mut status:Vec<u32> = vec![0;number];
    for i in 0..50000{
        //Save(&dna,&particles,&mut writer,i);
        particles = moveparticles(particles, &mut status, number);
        status = detectcollision(&particles, status, number, done, &dna, &mut boundto);
        boundpercent[i] = boundpercentfunction(&status, &number);
        if i%100 == 0{
            println!("{}",i);
            println!("{}/{}",boundpercent[i],number);
        }
    }
    saveboundpercent(boundpercent, &mut writer2);
    saveboundto(&mut boundto, &mut writer)?;
    Ok(())
}

fn placenewdna(dna: &mut Vec<Vec<f64>>, i:usize)->(Vec<f64>){
    let mut rng = rand::thread_rng();
    let mut newcords:Vec<f64> = vec![0.0,0.0];
    let mut end = false;
    while end == false{
        newcords[0] = rng.random::<f64>();
        newcords[1] = rng.random::<f64>();
        end = true;
        for i in 1..dna.len(){
            if (newcords[0]-dna[i][0]).powf(2.0) + (newcords[1]-dna[i][1]).powf(2.0) < RSQRD{
                end = false;
            }
        }
    }

    return newcords
}

fn placedna(dnanumber:usize)->(Vec<Vec<f64>>){
    let mut dna: Vec<Vec<f64>> = vec![vec![0.0;2];dnanumber];
    if dnanumber > 2{
        dna[0][0] = 0.5;
        dna[0][1] = 0.5;
    for i in 1..dnanumber{
        let newcords = placenewdna(&mut dna, i);
        dna[i][0] = newcords[0];
        dna[i][1] = newcords[1];
    }
    }
    else {
        dna[0][0] = 0.5;
        dna[0][1] = 0.5;
    }
    
    println!("{:?}",dna);
    return dna
}

fn Placeparticles(number:usize, dna: &mut Vec<Vec<f64>>)->Vec<Vec<f64>>{
    let mut particles: Vec<Vec<f64>> = vec![vec![0.0; 2]; number];
    let mut rng = rand::thread_rng();
    let centerx = 0.5;
    let centery = 0.5;
    let radius = 0.07;


    for i in 0..number{
        let angle = rng.gen::<f64>() * 2.0 * std::f64::consts::PI;
            
        particles[i][0] = centerx + radius * angle.cos();
        particles[i][1] = centery + radius * angle.sin();
        if i == 1{
            println!("particles[0] = {:?}", particles[0]);
            println!("particles[1] = {:?}", particles[1]);
        }
    }

    return particles
}

fn moveparticles(mut particles:Vec<Vec<f64>>, status:&mut Vec<u32>, number:usize)->Vec<Vec<f64>>{
    let mut rng = rand::thread_rng();
    for i in 0..number{
        if status[i] == 0 {
            particles[i][0] += rng.sample::<f64, _>(rand_distr::StandardNormal) / 100.0;
            particles[i][1] += rng.sample::<f64, _>(rand_distr::StandardNormal) / 100.0;
        if particles[i][0] > 1.0{
            particles[i][0] = 1.0;
        }
        if particles[i][0] < 0.0{
            particles[i][0] = 0.0;
        }
        if particles[i][1] > 1.0{
            particles[i][1] = 1.0;
        }
        if particles[i][1] < 0.0{
            particles[i][1] = 0.0;
        }
    }}


    return particles
}


fn detectcollision(particles:&Vec<Vec<f64>>, mut status:Vec<u32>, _number:usize, _done:usize, dna:&Vec<Vec<f64>>, boundto:&mut Vec<i32>)->Vec<u32>{
    let mut distsqrd:f64;
    
    for (particle_index, particle) in particles.iter().enumerate() {
        if status[particle_index] != 0 {
            continue; // Skip already bound particles
        }
        
        let mut smallest_dist = f64::MAX;
        let mut closest_strand = 0;
        let mut should_bind = false;
        
        // Find closest DNA strand
        for (strand_index, strand) in dna.iter().enumerate() {
            distsqrd = (particle[0]-strand[0]).powf(2.0)+(particle[1]-strand[1]).powf(2.0);
            if distsqrd < RSQRD && distsqrd < smallest_dist {
                smallest_dist = distsqrd;
                closest_strand = strand_index;
                should_bind = true;
            }
        }
        
        // Only try binding to closest strand
        if should_bind {
            let mut rng = rand::thread_rng();
            let random_number = rng.gen::<f64>();
            if random_number < BINDINGCHANCE { //Tal skal balanceres efter k_on
                status[particle_index] = closest_strand as u32 + 1; // Add 1 so 0 means unbound
                boundto[particle_index] = closest_strand as i32;
            }
        }
    }
    status
}

fn boundpercentfunction(status:&Vec<u32>, number:&usize)->i32{
    let mut amount:i32 = 0;
    for element in status{
        if *element != 0{
            amount+=1;
        }
    }
    return amount
}

fn testboundpercent(status:&Vec<u32>, number:&usize)->String{
    let mut amount:i32 = 0;
    for element in status{
        if *element != 0{
            amount+=1;
        }
    }
    let string_a:String = amount.to_string();
    let string_b:String = number.to_string();
    return string_a + " / " + &string_b;
}

fn Save(dna:&Vec<Vec<f64>>,particles:&Vec<Vec<f64>>, writer: &mut BufWriter<File>,i:usize)-> Result<(), Box<dyn std::error::Error>>{
    if i==0{ // First line of output.
        writeln!(writer,"DNA")?; //Write DNA as title. 
        for pair in dna.iter() {
            writeln!(writer, "{:.4},{:.4}", pair[0], pair[1])?; // Write "a,b" format for all DNA.
        }
        writeln!(writer, "Particles start");
        for tripple in particles.iter(){
            writeln!(writer, "{:.4},{:.4}",tripple[0],tripple[1])?; //Write the particle positions.
        }
    }
    else{
        writeln!(writer, "Particles {}",i)?;
        for tripple in particles.iter(){
            writeln!(writer, "{:.4},{:.4}",tripple[0],tripple[1])?; // Write the particle positions for every timestep after the first.
        }
    }
    
    writer.flush()?; //Flush the buffer.
    Ok(()) //Required by the function. 
}
fn saveboundpercent(boundpercent:Vec<i32>, writer2: &mut BufWriter<File>)-> Result<(), Box<dyn std::error::Error>>{
    for i in 0..boundpercent.len(){
        writeln!(writer2, "{}",boundpercent[i])?;
    }
    writer2.flush()?;
    Ok(())

}

fn saveboundto(boundto:&mut Vec<i32>, writer: &mut BufWriter<File>)-> Result<(), Box<dyn std::error::Error>>{
    for i in 0..boundto.len(){
        writeln!(writer, "{}",boundto[i])?;
    }
    writer.flush()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_moveparticles() {
        let mut particles = vec![vec![0.5, 0.5], vec![0.2, 0.2]];
        let mut status = vec![0, 0];
        let number = 2;

        let moved_particles = moveparticles(particles.clone(), &mut status, number);

        // Check that particles have moved
        assert_ne!(particles, moved_particles);
    }

    #[test]
    fn test_moveparticles_with_bound_particles() {
        let mut particles = vec![vec![0.5, 0.5], vec![0.2, 0.2]];
        let mut status = vec![1, 0]; // First particle is bound
        let number = 2;

        let moved_particles = moveparticles(particles.clone(), &mut status, number);

        // Check that bound particles have not moved
        assert_eq!(particles[0], moved_particles[0]);
        // Check that unbound particles have moved
        assert_ne!(particles[1], moved_particles[1]);
    }

    #[test]
    fn test_detectcollision() {
        let particles = vec![vec![0.5, 0.5], vec![0.2, 0.2]];
        let mut status = vec![0, 0];
        let number = 2;
        let done = 0;
        let dna = vec![vec![0.5, 0.5]];

        let updated_status = detectcollision(&particles, status.clone(), number, done, &dna, &mut vec![-1; number]);

        // Check that at least one particle has bound
        assert!(updated_status.iter().any(|&s| s != 0));
    }

    #[test]
    fn test_detectcollision_no_binding() {
        let particles = vec![vec![0.5, 0.5], vec![0.94, 0.94]];
        let mut status = vec![0, 0];
        let number = 2;
        let done = 0;
        let dna = vec![vec![1.0, 1.0]]; // DNA strand is far away

        let updated_status = detectcollision(&particles, status.clone(), number, done, &dna, &mut vec![-1; number]);

        // Check that no particles have bound
        assert!(updated_status.iter().all(|&s| s == 0));
    }

    #[test]
    fn test_detectcollision_multiple_dna_strands() {
        let particles = vec![vec![0.5, 0.5], vec![0.2, 0.2]];
        let mut status = vec![0, 0];
        let number = 2;
        let done = 0;
        let dna = vec![vec![0.5, 0.5], vec![0.2, 0.2]];
        let mut updated_status = detectcollision(&particles, status.clone(), number, done, &dna, &mut vec![-1; number]);
        for _ in 0..100{
            updated_status = detectcollision(&particles, status.clone(), number, done, &dna, &mut vec![-1; number]);
        }
        // Check that particles bind to the closest DNA strand
        assert_eq!(updated_status[0], 1); // First particle binds to first DNA strand
        assert_eq!(updated_status[1], 2); // Second particle binds to second DNA strand
    }

    #[test]
    fn test_boundpercent() {
        let status = vec![1, 0, 2, 0, 3];
        let number = 5;

        // Capture the output of boundpercent
        let output_str = testboundpercent(&status, &number);
        println!("restult = {}", output_str);
        assert_eq!(output_str, "3 / 5");
    }
}