use bevy::app::AppExit;
// use bevy::gizmos;
use bevy::{pbr::CascadeShadowConfigBuilder, prelude::*};
use bevy_editor_pls::prelude::*;
// use rand::Rng;
use std::f32::consts::PI;

mod latticebolt;
mod field;

use latticebolt::*;
use field::*;

fn main() {
    App::new()
        .add_state::<AppState>()
        .add_plugins(DefaultPlugins)
        .add_plugins(EditorPlugin::default())
        .add_plugins(FluidPlugin)
        .add_plugins(UnitPlugin)
        .run()
}

pub struct UnitPlugin;

impl Plugin for UnitPlugin {
    fn build(&self, app: &mut App) {
        app.add_systems(Startup, setup)
            .add_systems(Update, camera_update)
            .add_systems(Update, exit_game);
    }
}

pub struct FluidPlugin;

impl Plugin for FluidPlugin {
    fn build(&self, app: &mut App) {
        app.add_state::<FluidSimState>()
            .add_systems(
                FixedUpdate,
                (
                    apply_gravity,
                    // Needs multiple iterations per update
                    calc_divergence.after(apply_gravity),
                    solve_incompressibility.after(calc_divergence),
                ),
            )
            .add_systems(Update, sync_views.after(solve_incompressibility));
        app.insert_resource(Time::<Fixed>::from_seconds(1.0 / 60.0));
    }
}

#[derive(States, Debug, Clone, Copy, Eq, PartialEq, Hash, Default)]
pub enum FluidSimState {
    #[default]
    Reset,
    Running,
    Paused,
}

const GRAVITY: f32 = 9.81;
const COUNT: u32 = 10;
const SIZE: f32 = 2.0;
const CENTER: Vec3 = Vec3::new(
    (COUNT / 2) as f32 * SIZE,
    (COUNT / 2) as f32 * SIZE,
    (COUNT / 2) as f32 * SIZE,
);

const OVER_RELAX: f32 = 1.9;

#[derive(Component, Debug)]
struct MainCamera;

#[derive(Component, Debug)]
struct Velocity {
    pub v: Vec3,
    pub d: Vec3,
}

#[derive(Component, Debug)]
struct Divergence {
    pub v: Vec3,
    pub total: f32,
    pub c: f32,
}

#[derive(Component, Debug)]
struct Cell {
    pub pos: UVec3,
    pub open: f32,
}

#[derive(Resource, Debug)]
struct Config {
    pub size: UVec3,
    pub gravity: f32,
    pub grid_size: f32,
    pub over_relax: f32,
}

#[derive(States, Debug, Clone, Copy, Eq, PartialEq, Hash, Default)]
pub enum AppState {
    #[default]
    MainMenu,
    Game,
    GameOver,
}

fn setup(
    mut commands: Commands,
    mut meshes: ResMut<Assets<Mesh>>,
    mut materials: ResMut<Assets<StandardMaterial>>,
) {
    // commands.spawn(PbrBundle {
    //     mesh: meshes.add(shape::Circle::new(4.0).into()),
    //     material: materials.add(Color::WHITE.into()),
    //     transform: Transform::from_rotation(Quat::from_rotation_x(-std::f32::consts::FRAC_PI_2)),
    //     ..default()
    // });
    let config = Config {
        gravity: 9.81,
        size: UVec3::new(5, 5, 5),
        grid_size: 1.0,
        over_relax: 1.9,
    };

    commands.insert_resource(config);

    let mut _rng = rand::thread_rng();

    let mut field = Field::new(UVec3::new(COUNT, COUNT, 1));
    let center = field.size / 2;

    for x in 0..field.size.x {
        for y in 0..field.size.y {
            for z in 0..field.size.z {
                let cell = UVec3::new(x, y, z);
                let pos = Vec3::new(x as f32, y as f32, z as f32) * SIZE;
                let entity = commands
                    .spawn((
                        PbrBundle {
                            mesh: meshes.add(Mesh::from(shape::Box::new(0.1, 0.1, 1.0))),
                            material: materials.add(Color::rgba_u8(124, 144, 255, 128).into()),
                            transform: Transform::from_translation(pos),
                            ..default()
                        },
                        Velocity {
                            v: Vec3::ZERO,
                            d: Vec3::ZERO,
                            // v: Vec3::ZEROVec3::new(
                            //     rng.gen::<f32>() * 10.0,
                            //     rng.gen(),
                            //     rng.gen::<f32>() * 10.0,
                            // ),
                        },
                        Cell {
                            pos: cell,
                            open: if cell.y == 0
                                || cell.x == 0
                                || cell.y == field.size.y - 1
                                || cell.x == field.size.x - 1
                            {
                                0.0
                            } else {
                                1.0
                            },
                        },
                        Divergence {
                            v: Vec3::ZERO,
                            total: 0.0,
                            c: 0.0,
                        },
                    ))
                    .id();

                field.set_entity(cell, entity);
            }
        }
    }

    for (cell, entity) in field.cells.iter().enumerate() {
        println!("cell: {}, entity {}", cell, entity.index());
    }

    for z in 0..COUNT {
        for y in 0..COUNT {
            for x in 0..COUNT {
                println!(
                    "cell: ({}, {}, {}), entity: {}",
                    x,
                    y,
                    z,
                    field.get_entity(UVec3::new(x, y, z)).index()
                );
            }
        }
    }

    commands.insert_resource(field);

    commands.spawn(DirectionalLightBundle {
        directional_light: DirectionalLight {
            shadows_enabled: true,
            ..default()
        },
        transform: Transform {
            translation: Vec3::new(0.0, 2.0, 0.0),
            rotation: Quat::from_rotation_x(-PI / 4.),
            ..default()
        },

        cascade_shadow_config: CascadeShadowConfigBuilder {
            first_cascade_far_bound: 4.0,
            maximum_distance: 10.0,
            ..default()
        }
        .into(),
        ..default()
    });

    commands.spawn((
        Camera3dBundle {
            transform: Transform::from_xyz(-2.5, 4.5, 20.0).looking_at(CENTER, Vec3::Y),
            ..default()
        },
        MainCamera,
    ));

    commands.spawn(
        TextBundle::from_section(
            "A/D: Yaw, \nW/S: Pitch\nUp/Down: Distance",
            TextStyle {
                font_size: 20.0,
                ..default()
            },
        )
        .with_style(Style {
            position_type: PositionType::Absolute,
            top: Val::Px(12.0),
            left: Val::Px(12.0),
            ..default()
        }),
    );
}

fn camera_update(
    keyboard: Res<Input<KeyCode>>,
    mut q: Query<&mut Transform, With<MainCamera>>,
    time: Res<Time>,
) {
    let mut dir = Vec3::ZERO;

    if let Ok(mut transform) = q.get_single_mut() {
        if keyboard.pressed(KeyCode::A) {
            dir.y += -1.0;
        }

        if keyboard.pressed(KeyCode::D) {
            dir.y += 1.0;
        }

        if keyboard.pressed(KeyCode::Up) {
            dir.z += -1.0;
        }

        if keyboard.pressed(KeyCode::Down) {
            dir.z += 1.0;
        }

        if keyboard.pressed(KeyCode::W) {
            dir.x += -1.0;
        }

        if keyboard.pressed(KeyCode::S) {
            dir.x += 1.0;
        }

        if dir.length_squared() <= 0.0 {
            return;
        }

        dir *= time.delta_seconds();
        let forward = transform.local_z();
        transform.translation += forward * dir.z * 10.0;
        transform.rotate_around(
            CENTER,
            Quat::from_rotation_y(dir.y) * Quat::from_rotation_x(dir.x),
        );
        transform.look_at(CENTER, Vec3::Y);
    }
}

fn apply_gravity(mut q: Query<(&mut Velocity, &Cell)>, time: Res<Time<Fixed>>) {
    for (mut vel, cell) in q.iter_mut() {
        vel.v += Vec3::new(0.0, -GRAVITY, 0.0) * cell.open * time.delta_seconds();
    }
}

fn calc_divergence(
    mut q: Query<(&mut Divergence, &Cell, &Velocity)>,
    velocities: Query<&Velocity>,
    cells: Query<&Cell>,
    // time: Res<Time>,
    field: Res<Field>,
) {
    for (mut div, cell, vel) in q.iter_mut() {
        // var div = this.u[(i+1)*n + j] - this.u[i*n + j] +
        //		     this.v[i*n + j+1] - this.v[i*n + j];

        if cell.open == 0.0 {
            continue;
        }

        let mut d = Vec3::ZERO;
        let mut s: f32 = 0.0;

        if let Ok(rvel) = velocities.get(field.get_entity(cell.pos + UVec3::X)) {
            d.x += rvel.v.x - vel.v.x;
        }

        if let Ok(uvel) = velocities.get(field.get_entity(cell.pos + UVec3::Y)) {
            d.y += uvel.v.y - vel.v.y;
        }

        if let Ok(zvel) = velocities.get(field.get_entity(cell.pos + UVec3::Z)) {
            d.z += zvel.v.z - vel.v.z;
        }

        // Update solid status
        if let Ok(c) = cells.get(field.get_entity(cell.pos + UVec3::X)) {
            s += c.open;
        }

        if let Ok(c) = cells.get(field.get_entity(cell.pos + UVec3::Y)) {
            s += c.open;
        }

        if let Ok(c) = cells.get(field.get_entity(cell.pos + UVec3::Z)) {
            s += c.open;
        }

        if let Ok(c) = cells.get(field.get_entity(cell.pos.saturating_sub(UVec3::X))) {
            s += c.open;
        }

        if let Ok(c) = cells.get(field.get_entity(cell.pos.saturating_sub(UVec3::Y))) {
            s += c.open;
        }

        if let Ok(c) = cells.get(field.get_entity(cell.pos.saturating_sub(UVec3::Z))) {
            s += c.open;
        }

        if s == 0.0 {
            return;
        }
        // used to indicate which cells are solid

        div.v = d;
        div.total = d.x + d.y + d.z;
        div.c = -(div.total / s); // * OVER_RELAX;

        // println!("div{} = {}, vel{}", cell.pos, div.v, vel.v)
    }
}

fn solve_incompressibility(
    mut neighbors: Query<&Divergence>,
    mut q: Query<(&Divergence, &Cell, &mut Velocity)>,
    time: Res<Time<Fixed>>,
    field: Res<Field>,
) {
    for (div, cell, mut vel) in q.iter_mut() {
        if cell.open <= 0.0 {
            continue;
        }

        let mut d = Vec3::ZERO;
        let mut v = vel.v;
        // var div = this.u[(i+1)*n + j] - this.u[i*n + j] +
        //		     this.v[i*n + j+1] - this.v[i*n + j];
        if let Ok(rdiv) = neighbors.get_mut(field.get_entity(cell.pos.saturating_sub(UVec3::X))) {
            d.x += rdiv.c;
        }

        if let Ok(udiv) = neighbors.get_mut(field.get_entity(cell.pos.saturating_sub(UVec3::Y))) {
            d.y += udiv.c;
        }

        if let Ok(zdiv) = neighbors.get_mut(field.get_entity(cell.pos.saturating_sub(UVec3::Z))) {
            d.z += zdiv.c;
        }

        d *= time.delta_seconds();
        v -= div.c;
        v += d;
        vel.d = d;

        //vel.v = v.clamp_length(0.0, 100.0);

        //println!("{} c:{}, div{}, vel{}", cell.pos, div.c, div.v, vel.v)
    }
}

fn sync_views(mut q: Query<(&mut Transform, &Velocity, &Divergence, &Cell)>, mut gizmos: Gizmos) {
    for (mut transform, vel, _div, cell) in q.iter_mut() {
        transform.look_to(vel.v, vel.v.cross(Vec3::Y));
        if cell.open >= 1.0 {
            transform.scale.z = (vel.v.length_squared() / 100.0).abs().clamp(0.0001, 1.0);
        } else if cell.open < 1.0 {
            transform.scale *= 1.0 - cell.open.clamp(0.0, 1.0);
        }

        // gizmos.ray(
        //     transform.translation,
        //     vel.v,
        //     Color::Hsla {
        //         hue: (vel.v.length() / 1000.0) * 360.0,
        //         saturation: 1.0,
        //         lightness: 0.8,
        //         alpha: 1.0,
        //     },
        // );

        gizmos.ray(transform.translation, Vec3::X * vel.d.x * 100.0, Color::RED);
        gizmos.ray(transform.translation, Vec3::Y * vel.d.y * 100.0, Color::GREEN);
        gizmos.ray(transform.translation, Vec3::Z * vel.d.z * 100.0, Color::BLUE);
    }
}

fn exit_game(keyboard: Res<Input<KeyCode>>, mut app_exit: EventWriter<AppExit>) {
    if keyboard.just_pressed(KeyCode::Escape) {
        app_exit.send(AppExit);
    }
}
