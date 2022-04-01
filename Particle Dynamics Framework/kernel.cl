#pragma OPENCL EXTENSION cl_khr_global_int32_base_atomics : enable
#pragma OPENCL EXTENSION cl_khr_global_int32_extended_atomics : enable
#pragma OPENCL EXTENSION cl_khr_local_int32_base_atomics : enable
#pragma OPENCL EXTENSION cl_khr_local_int32_extended_atomics : enable

#define GETADRESS(t1,t2,size) ((t1) < (t2)) ? ((size * t1) + t2 - ((t1*(t1+1))/2)) : ((size * t2) + t1 - ((t2*(t2+1))/2))
#define GRID_TAIL_SIZE      32
#define BEAD_NHOOD_SIZE     256

typedef struct Vector3
{
    float x, y, z, w;
} my_Vector3;


typedef struct Bead
{
    my_Vector3    position;
    my_Vector3    image_position;
    my_Vector3    velocity;
    my_Vector3    force;
    my_Vector3    velocity_old;
    my_Vector3    force_old;    
} _Bead;            


typedef struct ParticleInfo
{        
    int          type;
    float        mass;
    float        mass_inv;
    float        lambda_dt_mass_inv;
    float        half_dt_mass_inv;
    float        half_dt2_mass_inv;

} _particle_info;



typedef struct GridCellNh
{
    int n[13];
} _gridCellNh;



typedef struct ParticleNh
{
    int n[BEAD_NHOOD_SIZE];
} _particleNh;



typedef struct InteractionNh
{
    float n[BEAD_NHOOD_SIZE];
} _interactionNh;


typedef struct GridCell
{
    int n[GRID_TAIL_SIZE];
} _GridCell;

typedef struct SimulationSettings
{
    int   typeNumber;
    int   units;
    float cutoff;
    float cutoff_inv;
    float dt;
    float half_dt;
    float half_dt2;
    float lambda_dt;
    float gamma;
    float sigma;
    int gridCellsX;
    int gridCellsY;
    int gridCellsZ;
    int gridCellsXY;
    float halfBoxSq;
    my_Vector3 gridStep;
    my_Vector3 boxSize;
} _settings;


typedef struct GridParameters
{
    int gridCellsX;
    int gridCellsY;
    int gridCellsZ;
    int gridCellsXY;
    float halfBoxSq;
    my_Vector3 gridStep;
    my_Vector3 boxSize;
} _parameters;



typedef struct Pair
{
    int type;
    int p1;
    int p2;
    float c1;
    float c2;
    float c3;
    float c4;
    float c5;
} _pair;


typedef struct SpringStructure                                                                    
{     
    int        type;    
    float      restLength;                                                                
    float      stiffness;                                                                
    float      damping;                                                                            
    int        p1;                                                            
    int        p2;    
                                                        
} _spring;                                                                        

typedef struct AngleBond
{
    int        type;
    int        b1;
    int        b2;
    int        b3;
    float      angle;
    float      c1;
    float      c2;
    float      c3;
} angle_bond;

inline void atomic_add_global(volatile __global float *addr, float val)
{
       union{
           unsigned int u32;
           float        f32;
       } next, expected, current;
       current.f32    = *addr;
       do{
          expected.f32 = current.f32;
           next.f32     = expected.f32 + val;
           current.u32  = atomic_cmpxchg( (volatile __global unsigned int *)addr, 
                               expected.u32, next.u32);
       } while( current.u32 != expected.u32 );
}



void kernel AdvancePosition(global _Bead* bead_buffer, global _particle_info* bead_info_buffer,  _settings settings, int beads_num)            
{                                                                        
    int bead_id = get_global_id(0);                                                                                                                        
    if(bead_id < beads_num)                                                            
    {    
        _Bead       bead      = bead_buffer[bead_id];
        float       lambda_dt = bead_info_buffer[bead_id].lambda_dt_mass_inv;
        float       half_dt2  = bead_info_buffer[bead_id].half_dt2_mass_inv;
        my_Vector3  increment = bead.velocity_old;

        bead.velocity.x = mad(bead.force.x, lambda_dt, increment.x);                    
        bead.velocity.y = mad(bead.force.y, lambda_dt, increment.y);                
        bead.velocity.z = mad(bead.force.z, lambda_dt, increment.z);

        increment.x *= settings.dt;
        increment.y *= settings.dt;
        increment.z *= settings.dt;

        increment.x = mad(bead.force.x, half_dt2, increment.x);
        increment.y = mad(bead.force.y, half_dt2, increment.y);
        increment.z = mad(bead.force.z, half_dt2, increment.z);

        bead.position.x += increment.x;
        bead.position.y += increment.y;
        bead.position.z += increment.z;

        bead.image_position.x += increment.x;
        bead.image_position.y += increment.y;
        bead.image_position.z += increment.z;
                     
        bead.force_old = bead.force;
        bead.force.x = 0.0f;
        bead.force.y = 0.0f;
        bead.force.z = 0.0f;
                 
        if(bead.position.z < 0.0f) bead.position.z = settings.boxSize.z - 0.000001f;
        if(bead.position.y < 0.0f) bead.position.y = settings.boxSize.y - 0.000001f;
        if(bead.position.x < 0.0f) bead.position.x = settings.boxSize.x - 0.000001f;
        
        if(bead.position.z > settings.boxSize.z) bead.position.z = 0.000001f;
        if(bead.position.y > settings.boxSize.y) bead.position.y = 0.000001f;
        if(bead.position.x > settings.boxSize.x) bead.position.x = 0.000001f;
        
        bead_buffer[bead_id] = bead;
    }                                                            
}                                                            


void kernel ResolveAngles(global angle_bond* bond, global _Bead* bead_buffer, int bond_num)
{
    int bond_id = get_global_id(0);                                                        
    if( bond_id < bond_num)                                                        
    {        
        angle_bond _bond = bond[bond_id];                                            
        
        my_Vector3 temp = bead_buffer[_bond.b2].image_position;                                                            
        
        my_Vector3 ba = bead_buffer[_bond.b1].image_position; 
        ba.x -= temp.x;
        ba.y -= temp.y;
        ba.z -= temp.z;
        
        my_Vector3 bc = bead_buffer[_bond.b3].image_position;                                               
        bc.x -= temp.x;
        bc.y -= temp.y;
        bc.z -= temp.z;
        
        float l_ba = sqrt(ba.x*ba.x + ba.y*ba.y + ba.z*ba.z);
        float l_bc = sqrt(bc.x*bc.x + bc.y*bc.y + bc.z*bc.z);
        
        if (l_ba != 0.0f) l_ba = 1.0f / l_ba;
        if (l_bc != 0.0f) l_bc = 1.0f / l_bc;
        
        ba.x *= l_ba;
        ba.y *= l_ba;
        ba.z *= l_ba;
        
        bc.x *= l_bc;
        bc.y *= l_bc;
        bc.z *= l_bc;
        
        float dotp = ba.x*bc.x + ba.y*bc.y + ba.z*bc.z;
        if(dotp < -1.0f) dotp = -0.9999f;
        if(dotp >  1.0f) dotp =  0.9999f;
        
        float angle = acos(dotp);
        float magnitude = -2.0f * _bond.c1 * (angle - _bond.angle);
        
        float fa = magnitude * l_ba;
        float fc = magnitude * l_bc;
        
        temp.x = ba.y*bc.z - ba.z*bc.y;
        temp.y = ba.z*bc.x - ba.x*bc.z;
        temp.z = ba.x*bc.y - ba.y*bc.x;
        
        //float ln = sqrt(temp.x*temp.x + temp.y*temp.y + temp.z*temp.z);
        //if (ln != 0.0f) ln = 1.0f / ln;
        //temp.x *= ln;
        //temp.y *= ln;
        //temp.z *= ln;
        
        //temp tez znormalizować
        my_Vector3 pa; 
        my_Vector3 pc; 
        
        pa.x = (ba.y*temp.z - ba.z*temp.y);
        pa.y = (ba.z*temp.x - ba.x*temp.z);
        pa.z = (ba.x*temp.y - ba.y*temp.x);
        
        float ln = sqrt(pa.x*pa.x + pa.y*pa.y + pa.z*pa.z);
        if (ln != 0.0f) ln = 1.0f / ln;
        ln *= fa;
        pa.x *= ln;
        pa.y *= ln;
        pa.z *= ln;
        
        
        pc.x = (temp.y*bc.z - temp.z*bc.y);
        pc.y = (temp.z*bc.x - temp.x*bc.z);
        pc.z = (temp.x*bc.y - temp.y*bc.x);
        
        ln = sqrt(pc.x*pc.x + pc.y*pc.y + pc.z*pc.z);
        if (ln != 0.0f) ln = 1.0f / ln;
        ln *= fc;
        pc.x *= ln;
        pc.y *= ln;
        pc.z *= ln;
        
        atomic_add_global(&(bead_buffer[_bond.b1].force.x), pa.x);
        atomic_add_global(&(bead_buffer[_bond.b1].force.y), pa.y);
        atomic_add_global(&(bead_buffer[_bond.b1].force.z), pa.z);
        
        atomic_add_global(&(bead_buffer[_bond.b3].force.x), pc.x);
        atomic_add_global(&(bead_buffer[_bond.b3].force.y), pc.y);
        atomic_add_global(&(bead_buffer[_bond.b3].force.z), pc.z);
        
        atomic_add_global(&(bead_buffer[_bond.b2].force.x), -pa.x-pc.x);
        atomic_add_global(&(bead_buffer[_bond.b2].force.y), -pa.y-pc.y);
        atomic_add_global(&(bead_buffer[_bond.b2].force.z), -pa.z-pc.z);
    }
}

void kernel ResolveSprings(global _Bead* bead_buffer, global _spring* springs, int springs_num)
{                                                                    
    int spring_id = get_global_id(0);                                                        
    if( spring_id < springs_num)                                                        
    {        
        _spring     spring  = springs[spring_id];                                                            
        int         p1      = spring.p1;                                                            
        int         p2      = spring.p2;                                                
        my_Vector3  ds      = bead_buffer[p1].image_position;                                                            
        my_Vector3  temp    = bead_buffer[p2].image_position;                                                                                                                    
        
        ds.x -= temp.x;                                
        ds.y -= temp.y;                                
        ds.z -= temp.z;                                

        float l      = native_sqrt(ds.x*ds.x + ds.y*ds.y + ds.z*ds.z);                        
        float factor = spring.restLength - l;                                
        if (l != 0.0f) l = 1.0f / l;                                        
        ds.x *= l;                                                
        ds.y *= l;                                                
        ds.z *= l;                                                
                       
        factor *= spring.stiffness;                                        
        ds.x*=factor;                                                
        ds.y*=factor;                                                
        ds.z*=factor;        
        
        atomic_add_global(&(bead_buffer[p1].force.x), ds.x);
        atomic_add_global(&(bead_buffer[p1].force.y), ds.y);
        atomic_add_global(&(bead_buffer[p1].force.z), ds.z);                                
        
        ds.x*=-1.0f;
        ds.y*=-1.0f;
        ds.z*=-1.0f;
        atomic_add_global(&(bead_buffer[p2].force.x), ds.x);
        atomic_add_global(&(bead_buffer[p2].force.y), ds.y);
        atomic_add_global(&(bead_buffer[p2].force.z), ds.z);
    }                                                    
}                                                    

        

void kernel AdvanceVelocity(global _Bead* bead_buffer, global _particle_info* bead_info_buffer, global _interactionNh* random_buffer, int beads_num, int buffer_size)
{                                                    
    int bead_id = get_global_id(0);                      
      
    _interactionNh random_value = random_buffer[bead_id];
    barrier(CLK_GLOBAL_MEM_FENCE);
     
    float temp = random_value.n[0];
    
    for(int ii = 0; ii < BEAD_NHOOD_SIZE-1; ii++)
    {
       random_value.n[ii] = random_value.n[ii+1];
    }
    random_value.n[BEAD_NHOOD_SIZE-1] = temp;
    random_buffer[bead_id] = random_value;
    
    barrier(CLK_GLOBAL_MEM_FENCE);
    
    if(bead_id < beads_num)                                        
    {    
        
        float      half_dt      = bead_info_buffer[bead_id].half_dt_mass_inv;
        my_Vector3 velocity_old = bead_buffer[bead_id].velocity_old;
        my_Vector3 force        = bead_buffer[bead_id].force;

        force.x += bead_buffer[bead_id].force_old.x;
        force.y += bead_buffer[bead_id].force_old.y;
        force.z += bead_buffer[bead_id].force_old.z;
        
        velocity_old.x = mad(force.x, half_dt, velocity_old.x);
        velocity_old.y = mad(force.y, half_dt, velocity_old.y);
        velocity_old.z = mad(force.z, half_dt, velocity_old.z);
        bead_buffer[bead_id].velocity_old = velocity_old;   
    }  
}



void kernel ResetBeadEnergy(global _Bead* bead_buffer, int beads_num)
{                                                    
    int bead_id = get_global_id(0);                                        
    if(bead_id < beads_num)                                        
    {    
        my_Vector3 empty_vector;
        empty_vector.x = 0.0f;
        empty_vector.y = 0.0f;
        empty_vector.z = 0.0f;

        bead_buffer[bead_id].force        = empty_vector;
        bead_buffer[bead_id].velocity     = empty_vector;
        bead_buffer[bead_id].force_old    = empty_vector;
        bead_buffer[bead_id].velocity_old = empty_vector;
    }
}



void kernel FillGrid(global _Bead* bead_buffer,  global int* grid_head, global _GridCell* grid_tail, _settings settings, int beads_num)
{                                                    
    int bead_id = get_global_id(0);                                        
    if(bead_id < beads_num)                                        
    {    
        my_Vector3 position = bead_buffer[bead_id].position;
        int cellX = floor((position.x) / settings.gridStep.x);
        int cellY = floor((position.y) / settings.gridStep.y);
        int cellZ = floor((position.z) / settings.gridStep.z);

        if( cellX >= 0 && cellX < settings.gridCellsX &&
            cellY >= 0 && cellY < settings.gridCellsY &&
            cellZ >= 0 && cellZ < settings.gridCellsZ)
        {
            int cell_id = cellX + cellY*settings.gridCellsX + cellZ*settings.gridCellsXY;
            //test if cell_id < grid_size
            int old = atomic_add(&(grid_head[cell_id]), 1);
            if(old < GRID_TAIL_SIZE) grid_tail[cell_id].n[old] = bead_id;
        }
    }
}

void kernel ResetGrid(global int* grid_head,  int grid_size)
{
    int grid_cell_id = get_global_id(0);
    if (grid_cell_id < grid_size)
    {
        grid_head[grid_cell_id] = 0;
    }
}


void kernel ResetContactList(global int* bead_nh_head, int bead_nh_list_size)
{
    int bead_id = get_global_id(0);
    if (bead_id < bead_nh_list_size)
    {
        bead_nh_head[bead_id] = 0;
    }
}


void kernel BuildContactList(global int* grid_head, global _GridCell* grid_tail, global _gridCellNh* grid_cell_nh, global int* bead_nh_head, global _particleNh* bead_nh_tail, global _Bead* bead_buffer, _settings settings, int grid_size)
{
    int grid_cell_id = get_global_id(0);
    if (grid_cell_id < grid_size)
    {
        int beads_in_cell = grid_head[grid_cell_id];
        if( beads_in_cell > 0 )
        {
            if(beads_in_cell > GRID_TAIL_SIZE) beads_in_cell = GRID_TAIL_SIZE;
            
            _GridCell    bead_id = grid_tail[grid_cell_id];
            _GridCell    nh_bead_id;
            
            for(int ii = 0; ii < (beads_in_cell-1); ii++)
            {
                for(int jj = (ii+1); jj < beads_in_cell; jj++)
                {
                    int size = atomic_add(&(bead_nh_head[bead_id.n[ii]]), 1);
                    if( size < BEAD_NHOOD_SIZE) bead_nh_tail[bead_id.n[ii]].n[size] = bead_id.n[jj];
                }    
            }
            
            _gridCellNh nh = grid_cell_nh[grid_cell_id];
            for(int ii = 0; ii < 13; ii++)    
            {
                int nh_cell_id = nh.n[ii];
                if( nh_cell_id != -1 )
                {
                    int nh_beads_in_cell = grid_head[nh_cell_id];
                    if( nh_beads_in_cell > 0 )
                    {
                        nh_bead_id = grid_tail[nh_cell_id];
                        
                        if( nh_beads_in_cell > GRID_TAIL_SIZE ) nh_beads_in_cell = GRID_TAIL_SIZE;
                        for(int jj = 0; jj < nh_beads_in_cell; jj++)
                        {
                            for(int kk = 0; kk < (beads_in_cell); kk++)
                            {
                                int size = atomic_add(&(bead_nh_head[bead_id.n[kk]]), 1);
                                if( size < BEAD_NHOOD_SIZE) bead_nh_tail[bead_id.n[kk]].n[size] = nh_bead_id.n[jj];
                            }
                        }
                    }
                }
            }     
        }
    }
}



void kernel InitContactList(global int* bead_nh_head, global _particleNh* bead_nh_tail, global _particle_info* bead_info, global _particleNh* interaction_id , __global const _pair* coefficient_table, int interaction_types_num, int beads_num)
{
    int bead_id = get_global_id(0);
    if (bead_id < beads_num)
    {
        int number_of_neigbours = bead_nh_head[bead_id];
        if(number_of_neigbours > 0)
        {
            if(number_of_neigbours > BEAD_NHOOD_SIZE) number_of_neigbours = BEAD_NHOOD_SIZE;
            int type_1 = bead_info[bead_id].type;

            for(int ii = 0; ii < number_of_neigbours; ii++)
            {    
                int type_2  = bead_info[bead_nh_tail[bead_id].n[ii]].type;
                int pair_id = GETADRESS(type_1,type_2,interaction_types_num);
                
                interaction_id[bead_id].n[ii] = pair_id;
            }
        }
    }
}


void kernel ResolveContacts(global int* bead_nh_head, global _particleNh* bead_nh_tail, global _Bead* bead_buffer, global _interactionNh* random_buffer, global _particleNh* interaction_id, __global const _pair* coefficient_table, _settings settings, int beads_num)
{
    int bead_id_1 = get_global_id(0);
    
    if (bead_id_1 < beads_num)
    {
        int number_of_neigbours = bead_nh_head[bead_id_1];
        if(number_of_neigbours > 0)
        {
            if(number_of_neigbours > BEAD_NHOOD_SIZE) number_of_neigbours = BEAD_NHOOD_SIZE;
            
            float cutoff_sq = 1.0f;
            if (settings.units != 0)
            {
                 cutoff_sq = settings.cutoff * settings.cutoff;
            }
            
            int   pair_id   = -1;
            float alpha     = 0.0f;
            float rd        = 0.0f;
            float gamma     = settings.gamma;
            float sigma     = settings.sigma;
            float maxDim    = settings.halfBoxSq;
            float random    = 0.0f;
               
            my_Vector3 Reaction;
            Reaction.x = 0.0f;
            Reaction.y = 0.0f;
            Reaction.z = 0.0f;

            my_Vector3 boxSize    = settings.boxSize;
            my_Vector3 position_1 = bead_buffer[bead_id_1].position;
            my_Vector3 velocity_1 = bead_buffer[bead_id_1].velocity;

            my_Vector3 position_2;
            my_Vector3 velocity_2;
            my_Vector3 rij;
            my_Vector3 vij;
        
            for(int ii = 0; ii < number_of_neigbours; ii++)
            {   
                int bead_id_2 = bead_nh_tail[bead_id_1].n[ii];
                position_2 = bead_buffer[bead_id_2].position;

                rij.x = position_1.x - position_2.x;
                rij.y = position_1.y - position_2.y;
                rij.z = position_1.z - position_2.z;
    
                float length  = rij.x*rij.x + rij.y*rij.y + rij.z*rij.z;

                if(length > maxDim)
                {
                    if(rij.x*rij.x > maxDim)
                    {
                        if(rij.x > 0.0f)
                            rij.x -= boxSize.x;
                        else    
                            rij.x += boxSize.x;
                     }

                    if(rij.y*rij.y > maxDim)
                    {
                        if(rij.y > 0.0f)
                            rij.y -= boxSize.y;
                        else    
                            rij.y += boxSize.y;
                    }

                    if(rij.z*rij.z > maxDim)
                    {
                        if(rij.z > 0.0f)
                            rij.z -= boxSize.z;
                        else    
                            rij.z += boxSize.z;
                    }
                    length = rij.x*rij.x + rij.y*rij.y + rij.z*rij.z;
                }


                if(length < cutoff_sq) 
                {        
                    pair_id = interaction_id[bead_id_1].n[ii];
                    alpha = coefficient_table[pair_id].c1;
                    rd    = coefficient_table[pair_id].c2;
                    length = sqrt(length);
                    
                    float coeff = rd*length + 1.0f;
                    
                    if (coeff > 0.0f)
                    {
                        random = random_buffer[bead_id_1].n[ii]; 
                        

                        if(settings.units != 0) length *= settings.cutoff_inv;

                        if(length != 0.0f)
                            length = 1.0f/length;
                        else
                            length = 0.0f;

                        rij.x *= length;
                        rij.y *= length;
                        rij.z *= length;

                        velocity_2 = bead_buffer[bead_id_2].velocity;
                        vij.x = velocity_1.x - velocity_2.x;
                        vij.y = velocity_1.y - velocity_2.y;
                        vij.z = velocity_1.z - velocity_2.z;
                
                        float scalar = vij.x*rij.x + vij.y*rij.y + vij.z*rij.z;

                        float fij =   alpha * coeff;
                        //fij -= gamma * coeff * coeff * scalar;
                        //fij += sigma * random * coeff;
                
                        float wd = pow(coeff,0.5f);
                        fij -= gamma * wd * scalar;
                        fij += sigma * random * sqrt(wd);

                        Reaction.x = fij * rij.x;
                        Reaction.y = fij * rij.y;
                        Reaction.z = fij * rij.z;

                        atomic_add_global(&(bead_buffer[bead_id_1].force.x), Reaction.x);
                        atomic_add_global(&(bead_buffer[bead_id_1].force.y), Reaction.y);
                        atomic_add_global(&(bead_buffer[bead_id_1].force.z), Reaction.z);

                        atomic_add_global(&(bead_buffer[bead_id_2].force.x), -Reaction.x);
                        atomic_add_global(&(bead_buffer[bead_id_2].force.y), -Reaction.y);
                        atomic_add_global(&(bead_buffer[bead_id_2].force.z), -Reaction.z);
                    }
                }
            }
        }
    }
  
}






