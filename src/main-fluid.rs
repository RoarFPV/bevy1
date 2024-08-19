use bevy::app::AppExit;
// use bevy::gizmos;
use bevy::{pbr::CascadeShadowConfigBuilder, prelude::*};
use bevy_editor_pls::prelude::*;
use latticebolt::LatticeConfig;
use rand::Rng;
// use rand::Rng;
use std::f32::consts::PI;

mod field;
mod latticebolt;

use field::*;

fn main() {
    App::new()
        .add_plugins(DefaultPlugins)
        // .add_plugins(EditorPlugin::default())
        .add_state::<AppState>()
        .add_state::<FluidSimState>()
        .insert_resource(Time::<Fixed>::from_seconds(1.0 / 60.0))
        .insert_resource(Field::new(IVec3::new(50, 20, 1)))
        .insert_resource(LatticeConfig::default())
        .add_systems(Startup, setup)
        .add_systems(
            FixedUpdate,
            (latticebolt::stream, latticebolt::collide)
                .chain()
                .run_if(in_state(FluidSimState::Running)),
        )
        .add_systems(Update, camera_update)
        .add_systems(Update, exit_game)
        .add_systems(Update, sync_views.after(latticebolt::collide))
        .run()
}

#[derive(States, Debug, Clone, Copy, Eq, PartialEq, Hash, Default)]
pub enum FluidSimState {
    #[default]
    Reset,
    Running,
    Paused,
}

#[derive(Component, Debug)]
struct MainCamera;

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
    mut field: ResMut<Field>,
    lattice_config: Res<LatticeConfig>,
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

    let mut rng = rand::thread_rng();

    for x in 0..field.size.x {
        for y in 0..field.size.y {
            for z in 0..field.size.z {
                let cell = IVec3::new(x, y, z);
                let pos = Vec3::new(x as f32, y as f32, z as f32) * config.grid_size;
                let entity = commands
                    .spawn((
                        // PbrBundle {
                        //     mesh: meshes.add(Mesh::from(shape::Box::new(0.1, 0.1, 1.0))),
                        //     material: materials.add(Color::rgba_u8(124, 144, 255, 128).into()),
                        //     transform: Transform::from_translation(pos),
                        //     ..default()
                        // },
                        latticebolt::Cell {
                            pos: cell,
                            velocity: Vec3::ZERO,
                            density: [
                                if cell.x == field.size.x / 2 + 1 && cell.y == field.size.y / 2 {
                                    lattice_config.rho0 * 10.0
                                } else {
                                    lattice_config.rho0 + 0.01 * rng.gen::<f32>()
                                },
                                lattice_config.rho0 + 0.01 * rng.gen::<f32>(),
                                lattice_config.rho0 + 0.01 * rng.gen::<f32>(),
                                lattice_config.rho0 + 0.01 * rng.gen::<f32>(),
                                lattice_config.rho0 + 0.01 * rng.gen::<f32>(),
                                lattice_config.rho0 + 0.01 * rng.gen::<f32>(),
                                lattice_config.rho0 + 0.01 * rng.gen::<f32>(),
                                lattice_config.rho0 + 0.01 * rng.gen::<f32>(),
                                lattice_config.rho0 + 0.01 * rng.gen::<f32>(),
                            ],
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
            transform: Transform::from_xyz(-2.5, 4.5, 20.0).looking_at(field.center(), Vec3::Y),
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

    commands.insert_resource(config);
}

fn camera_update(
    keyboard: Res<Input<KeyCode>>,
    mut q: Query<&mut Transform, With<MainCamera>>,
    field: Res<Field>,
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
        let center = field.size.as_vec3() / 2.0;
        transform.translation += forward * dir.z * 10.0;
        transform.rotate_around(
            center,
            Quat::from_rotation_y(dir.y) * Quat::from_rotation_x(dir.x),
        );
        transform.look_at(center, Vec3::Y);
    }
}

fn sync_views(mut q: Query<&latticebolt::Cell>, mut gizmos: Gizmos) {
    for cell in q.iter_mut() {
        let vel = cell.velocity;
        //transform.look_to(vel, vel.cross(Vec3::Y));

        //transform.scale.z = (cell.density[0] / 100.0).clamp(0.0, 2.0);
        // if cell.density[0] >= 1.0 {
        //     transform.scale.z = (vel.length_squared() / 100.0).abs().clamp(0.0001, 1.0);
        // } else if cell.density[0] < 1.0 {
        //     transform.scale *= 1.0 - cell.density[0].clamp(0.0, 1.0);
        // }

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

        // gizmos.ray(transform.translation, vel, Color::RED);

        gizmos.ray(
            cell.pos.as_vec3(),
            vel.normalize() / 2.0,
            Color::hsla((cell.density[0] / 5.0) * 360.0, 0.7, 0.7, 0.5),
        );
        // gizmos.ray(transform.translation, Vec3::Y * vel.d.y * 100.0, Color::GREEN);
        // gizmos.ray(transform.translation, Vec3::Z * vel.d.z * 100.0, Color::BLUE);
    }
}

fn exit_game(keyboard: Res<Input<KeyCode>>, mut app_exit: EventWriter<AppExit>) {
    if keyboard.just_pressed(KeyCode::Escape) {
        app_exit.send(AppExit);
    }
}
