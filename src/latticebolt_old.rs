use std::{sync::Arc};

use bevy::{
    ecs::query::{BatchingStrategy, QueryEntityError},
    math::{IVec3, UVec3, Vec3},
    //utils::hashbrown::{hash_map, HashMap},
    prelude::*,
    tasks::ComputeTaskPool,
    utils::{petgraph::visit::IntoNeighbors, Hashed, PreHashMap},
};

pub const NODE_COUNT: usize = 19;
pub struct LaticeBoltPlugin;

impl Plugin for LaticeBoltPlugin {
    fn build(&self, app: &mut App) {
        let config = LatticeConfig::default();

        app.add_systems(FixedUpdate, (step).chain());
        app.insert_resource(config);
    }
}

#[derive(Component, Debug, Copy, Clone, Reflect)]
pub struct Cell {
    pub id: CellIndex,
    pub pos: IVec3,
    pub density: [f32; NODE_COUNT],
    pub total_density: f32,
    pub velocity: Vec3,
}

#[derive(Debug, Default, Copy, Clone, Reflect)]
pub struct Node {
    pub weight: f32,
    pub velocity: Vec3,
}

impl Node {
    fn new(weight: f32, velocity: Vec3) -> Self {
        Self { weight, velocity }
    }
}

#[derive(Resource, Debug, Reflect, Copy, Clone)]
#[reflect(Resource)]
pub struct LatticeConfig {
    // average density
    pub rho0: f32,
    // collision timescale
    pub tau: f32,

    

    pub density_max: f32,

    pub view_rho: f32,
    pub view_density_max: f32,
    pub eq: f32,
    pub size: UVec3,
    pub emitter: UVec3,

    // lattice velocity and weights
    pub nodes: [Node; 19],

    
    
}

const REFLECTED_NODES: [usize; NODE_COUNT] = [
    0,  // 0 -> 0, center -> center
    3,  // 1 -> 3, north -> south
    4,  // 2 -> 4, east -> west
    1,  // 3 -> 1, south -> north
    2,  // 4 -> 2, west -> east
    6,  // 5 -> 6, top -> bottom
    5,  // 6 -> 5, bottom -> top
    9,  // 7 -> 9, north east -> south west
    10, // 8 -> 10, south east -> north west
    7,  // 9 -> 7, south west -> north east
    8,  // 10 -> 8, north west -> south east
    17, // 11 -> 17, top north -> bottom south
    18, // 12 -> 18, top east -> bottom west
    15, // 13 -> 15, top south -> bottom north
    16, // 14 -> 16, top west -> bottom east
    13, // 15 -> 13, bottom north -> top south
    14, // 16 -> 14, bottom east -> top west
    11, // 17 -> 11, bottom south -> top north
    12, // 18 -> 12, bottom west -> top east
];

impl Default for LatticeConfig {
    fn default() -> Self {
        Self {
            rho0: 100.0,
            view_rho: 0.0001,
            view_density_max: 10.0,
            tau: 1.99,
            density_max: 1000.0,
            eq: 1.0,
            nodes: [
                Node::new(1.0 / 3.0, Vec3::new(0.0, 0.0, 0.0)),
                Node::new(1.0 / 18.0, Vec3::new(0.0, 1.0, 0.0)), // north
                Node::new(1.0 / 18.0, Vec3::new(1.0, 0.0, 0.0)), // east
                Node::new(1.0 / 18.0, Vec3::new(0.0, -1.0, 0.0)), // south
                Node::new(1.0 / 18.0, Vec3::new(-1.0, 0.0, 0.0)), // west
                Node::new(1.0 / 18.0, Vec3::new(0.0, 0.0, 1.0)), // top
                Node::new(1.0 / 18.0, Vec3::new(0.0, 0.0, -1.0)), // bottom
                Node::new(1.0 / 36.0, Vec3::new(1.0, 1.0, 0.0)), // north east
                Node::new(1.0 / 36.0, Vec3::new(1.0, -1.0, 0.0)), // south east
                Node::new(1.0 / 36.0, Vec3::new(-1.0, -1.0, 0.0)), // south west
                Node::new(1.0 / 36.0, Vec3::new(-1.0, 1.0, 0.0)), // north west
                Node::new(1.0 / 36.0, Vec3::new(0.0, 1.0, 1.0)), // top north
                Node::new(1.0 / 36.0, Vec3::new(1.0, 0.0, 1.0)), // top east
                Node::new(1.0 / 36.0, Vec3::new(0.0, -1.0, 1.0)), // top south
                Node::new(1.0 / 36.0, Vec3::new(-1.0, 0.0, 1.0)), // top west
                Node::new(1.0 / 36.0, Vec3::new(0.0, 1.0, -1.0)), // bottom north
                Node::new(1.0 / 36.0, Vec3::new(1.0, 0.0, -1.0)), // bottom east
                Node::new(1.0 / 36.0, Vec3::new(0.0, -1.0, -1.0)), // bottom south
                Node::new(1.0 / 36.0, Vec3::new(-1.0, 0.0, -1.0)), // bottom west
            ],
            size: UVec3::new(5, 5, 1),
            emitter: UVec3::ZERO,
        }
    }
}

#[derive(Default, Resource)]
pub struct LatticeRuntime {
    pub gen: bool,
    pub id_to_entity: PreHashMap<CellIndex, Entity>,
}

pub type CellMap = PreHashMap<CellIndex, Cell>;
pub type CellIndex = u64;

#[derive(Default, Resource)]
pub struct LatticeRuntimeCurrent {
    pub cells: CellMap,
}

#[derive(Default, Resource)]
pub struct LatticeRuntimeNext {
    pub cells: CellMap,
}

impl LatticeRuntime {
    pub fn index(size: IVec3, pos: IVec3) -> CellIndex {
        // if pos.cmpge(size).any() || pos.cmplt(IVec3::ZERO).any() {
        //     return usize::MAX;
        // }

        if pos.x < 0
            || pos.x >= size.x
            || pos.y < 0
            || pos.y >= size.y
            || pos.z < 0
            || pos.z >= size.z
        {
            return CellIndex::MAX;
        }

        (pos.x + (size.x * pos.y) + (size.x * size.y * pos.z)) as CellIndex
    }

    pub fn next_gen(&mut self) {
        self.gen = !self.gen;
    }
}

pub fn step(
    config: Res<LatticeConfig>,
    mut lattice: ResMut<LatticeRuntime>,
    time: Res<Time<Fixed>>,
    mut res_cells: ResMut<LatticeRuntimeCurrent>,
    mut res_cells_next: ResMut<LatticeRuntimeNext>,
) {
    // start processing
    // make a look up table for neighbors
    let cells_next: &mut CellMap;
    let cells: &CellMap;

    if lattice.gen {
        cells_next = &mut res_cells.cells;
        cells = &res_cells_next.cells;
    } else {
        cells_next = &mut res_cells_next.cells;
        cells = &res_cells.cells;
    };

    let pool = ComputeTaskPool::get();

    let config_ = *config;
    let size = config.size.as_ivec3();

    let BATCH_COUNT: usize = (cells.len() * 2) / 200;
    let batch = cells.len() / BATCH_COUNT;

    let cells_arc = Arc::new(cells);

    let output = pool.scope(move |scope| {
        // creates a span and starts the timer
        for b in 0..batch + 1 {
            let start = b * BATCH_COUNT;

            let cells_scope = cells_arc.clone();
            scope.spawn(async move {
                // let _span = info_span!("batch {}", b, name = "batch_").entered();
                let mut nextcells = Vec::<Cell>::with_capacity(BATCH_COUNT);
                for bt in start..(start + BATCH_COUNT) {
                    if let Some((id, rcell)) = cells_scope.iter().nth(bt) {
                        nextcells.push(apply_collision_old(
                            rcell,
                            cells,
                            &config_,
                            1.0 / config_.tau,
                        ));
                    }
                }
                nextcells
            });
        }
    });

    {
        // let _span = info_span!("flatten", name = "lattice_flatten").entered();

        for cell in output.iter().flatten() {
            cells_next.insert(Hashed::new(cell.id), *cell);
        }
    }

    lattice.next_gen();
}

// // https://github.com/kynan/firesim/blob/7241cbfd34e10dd80ba78696c84ff32f976893e7/src/lbm/LBM_def.h#L892
// fn apply_collision_fsim(
//     cell: &Cell,
//     index:usize,
//     cells_next: &mut Vec<Cell>,
//     config: &LatticeConfig,
//     dt:f32,
// ) {
//     let mut rho = cell.density[0];
//     let mut total_vel = Vec3::ZERO;

//     for d in 1..config.nodes.len() {
//         let fi = cell.density[d];
//         rho += fi;
//         total_vel += config.nodes[d].velocity * fi;
//     }

//     // collision
//     let fc = rho - 1.5 * (total_vel.length_squared());
//     let omegai = 1.0 - config.tau;

//     cells_next[index].density[0] = omegai * cell.density[0] + config.tau * config.nodes[0].weight * fc;
//     cells_next[index].velocity = total_vel;

//     for d in 1..config.nodes.len() {
//         let node = &config.nodes[d];
//         let eiu = node.velocity.dot(total_vel);
//         let density = cell.density[d];

//         let next_index = LatticeRuntime::index(config.size.as_ivec3(), cell.pos + node.velocity.as_ivec3());

//         let next_density=omegai * density
//         + config.tau * node.weight * (fc + 3.0 * eiu + 4.5 * eiu * eiu);

//         if let Some(stream_cell) = cells_next.get_mut(next_index) {
//             stream_cell.density[d] = next_density;

//         } else {
//             cells_next[index].density[d] -= next_density;
//          }

//     }
// }

fn apply_collision_old(cell: &Cell, cells: &CellMap, config: &LatticeConfig, dt: f32) -> Cell {
    // let _span = info_span!("apply_collision", cell.id).entered();

    let mut next = *cell;

    for i in 1..cell.density.len() {
        let node = &config.nodes[i];
        // stream from outside
        let prev_index =
            LatticeRuntime::index(config.size.as_ivec3(), cell.pos - node.velocity.as_ivec3());
        if let Some(stream_cell) = cells.get(&Hashed::new(prev_index)) {
            let density = next.density[i];

            //next.density[i] += (stream_cell.density[i] - next.density[i])*dt;
            next.density[i] += (stream_cell.density[i]) * dt;

            // make sure to account for density streamed out in this direction
            next.density[i] -= (density) * dt;
        } else {
            let r = REFLECTED_NODES[i];
            let density = next.density[i];

            let delta = (next.density[r]) * dt;
            next.density[i] += delta;
            next.density[i] -= density * dt;
        }
    }

    let mut total_density: f32 = 0.0;
    for i in 0..next.density.len() {
        total_density += next.density[i]
    }

    next.total_density = total_density;

    let mut vel = Vec3::ZERO;
    for (i, node) in config.nodes.iter().enumerate() {
        vel += node.velocity * next.density[i];
    }

    //println!("vel: {}, d: {:?}", vel, next.density);
    vel /= total_density;

    // add any outside accelerations/velocities here
    //

    // gravity
    //vel += -Vec3::Y * 9.81 * dt;

    // collision
    let velsq = vel.length_squared();
    let mut eqlib = [0.0; NODE_COUNT];

    next.velocity = vel;

    // collision function
    for d in 0..config.nodes.len() {
        let n = config.nodes[d];
        let term = n.velocity.dot(vel);
        eqlib[d] = n.weight
            * total_density
            * (1.0 + 3.0 * term + 9.0 * (term * term) / 2.0 + -3.0 * velsq / 2.0)
            * config.rho0;

        let density = next.density[d];
        next.density[d] -= (density - eqlib[d]) * dt;
        next.density[d] = next.density[d].clamp(-config.density_max, config.density_max);
    }

    next
}

#[derive(Default, Copy, Clone)]
pub struct CellDensity {
    pub density: [f32; NODE_COUNT],
    pub total_density: f32,
    pub total_velocity: Vec3,
}

#[derive(Component, Default)]
// #[reflect(Component)]
pub struct CellDensityGenA {
    pub v: CellDensity,
}

#[derive(Component, Default)]
// #[reflect(Component)]
pub struct CellDensityGenB {
    pub v: CellDensity,
}

#[derive(Component, Default, Reflect, Copy, Clone)]
#[reflect(Component)]
pub struct CellId {
    pub id: CellIndex,
    pub pos: IVec3,
}

pub fn step_cell<F>(
    dt: f32,
    id_to_entity: &PreHashMap<CellIndex, Entity>,
    config: &LatticeConfig,
    read_id: &CellId,
    read: &CellDensity,
    write: &mut CellDensity,
    neighbors: F,
) where
    F: Fn(&Entity) -> Result<(CellId, CellDensity), QueryEntityError>
{
    for i in 1..read.density.len() {
        let node = &config.nodes[i];
        // stream from outside
        let prev_index =
            LatticeRuntime::index(config.size.as_ivec3(), read_id.pos - node.velocity.as_ivec3());

        write.density[i] = read.density[i];

        if let Some(neighbor_ent) = id_to_entity.get(&Hashed::new(prev_index)) {
            if let Ok((n_id, n_density)) = neighbors(neighbor_ent) {
                let density = read.density[i];

                write.density[i] = n_density.density[i];
                //next.density[i] += (stream_cell.density[i] - next.density[i])*dt;
                // write.density[i] += (n_density.density[i]) * dt;

                // // make sure to account for density streamed out in this direction
                // write.density[i] -= (density) * dt;
                continue;
            }
        }
        // not found so reflect
        let r = REFLECTED_NODES[i];
        write.density[i] = read.density[r];
        // let density = read.density[i];

        // let delta = (write.density[r]) * dt;
        // write.density[i] += delta;
        // write.density[i] -= density * dt;
    }

    let mut total_density: f32 = 0.0;
    for i in 0..read.density.len() {
        total_density += read.density[i]
    }

    //total_density = total_density.clamp(-config.density_max, config.density_max);

    write.total_density = total_density;

    let mut vel = Vec3::ZERO;
    for (i, node) in config.nodes.iter().enumerate() {
        vel += node.velocity * write.density[i];
    }

    //println!("vel: {}, d: {:?}", vel, next.density);
    vel /= total_density;

    vel.clamp_length_max(config.density_max);

    // add any outside accelerations/velocities here
    //

    // heat rising
    //vel += Vec3::Y * (write.density[0].abs()/(config.density_max*config.density_max)) ;

    // collision
    let velsq = vel.length_squared();
    

    write.total_velocity = vel;

    let h = 0.5+ (read_id.pos.y as f32/ config.size.y as f32); 
    // collision function
    for d in 0..config.nodes.len() {
        let n = config.nodes[d];
        let term = n.velocity.dot(vel);
        let eqlib = n.weight
            * total_density
            * (1.0 
                + (3.0 * term )
                + (9.0/2.0) * term * term
                + (-3.0/2.0) * velsq);
   
        let density = write.density[d];
        //write.density[d] -= (density - eqlib) * dt;

        let f = (1.0 - dt) * density + dt * eqlib;

        write.density[d] = f;

        //write.density[d] = write.density[d].clamp(-config.density_max, config.density_max);
    }
}

pub fn step_gen_a_to_b(
    config: Res<LatticeConfig>,
    mut runtime: ResMut<LatticeRuntime>,
    mut read: Query<(&CellId, &CellDensityGenA, &mut CellDensityGenB)>,
    neighbors: Query<(&CellId, &CellDensityGenA)>,
    time: Res<Time<Fixed>>,
) {
    let gen = runtime.gen;
    if !gen {
        return;
    }

    let dt = 1.0/config.tau;
    let nr = Arc::new(neighbors);
    read
        //.iter_mut()
        .par_iter_mut()
        .batching_strategy(BatchingStrategy::fixed(128))
        .for_each(|(id, read_gen, mut write_gen)| {
            let n = nr.clone();
            step_cell(
                dt, 
                &runtime.id_to_entity, 
                &config,
                id, 
                &read_gen.v, 
                &mut write_gen.v,
                 move |e| {
                if let Ok((id, cellgen)) = n.get(*e) {
                    Ok((*id, cellgen.v))
                } else {
                    Err(QueryEntityError::NoSuchEntity(*e))
                }
            });
        });
}


pub fn step_gen_b_to_a(
    config: Res<LatticeConfig>,
    mut runtime: ResMut<LatticeRuntime>,
    mut read: Query<(&CellId, &CellDensityGenB, &mut CellDensityGenA)>,
    neighbors: Query<(&CellId, &CellDensityGenB)>,
    time: Res<Time<Fixed>>,
) {
    let gen = runtime.gen;
    if gen {
        return;
    }


    let dt = 1.0/config.tau;
    let nr = Arc::new(neighbors);
    read
        //.iter_mut()
        .par_iter_mut()
        .batching_strategy(BatchingStrategy::fixed(128))
        .for_each(|(id, read_gen, mut write_gen)| {
            let n = nr.clone();
            step_cell(
                dt, 
                &runtime.id_to_entity, 
                &config,
                id, 
                &read_gen.v, 
                &mut write_gen.v,
                 move |e| {
                if let Ok((id, cellgen)) = n.get(*e) {
                    Ok((*id, cellgen.v))
                } else {
                    Err(QueryEntityError::NoSuchEntity(*e))
                }
            });
        });
}


pub fn next_gen(mut runtime: ResMut<LatticeRuntime>) {
    let gen = runtime.gen;
    runtime.gen = !gen;
}