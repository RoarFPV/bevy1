use bevy::app::AppExit;
use bevy::prelude::*;

use crate::ui::State;
use crate::ui::HOVERED_BUTTON;
use crate::ui::NORMAL_BUTTON;
use crate::ui::PRESSED_BUTTON;

#[derive(Debug)]
enum ButtonActions {
    Play,
    Quit,
}

#[derive(Component, Debug)]
struct ButtonAction {
    action: ButtonActions,
}

#[derive(Resource)]
struct MenuData {
    root: Entity,
}

pub struct MenuPlugin;

impl Plugin for MenuPlugin {
    fn build(&self, app: &mut App) {
        app.add_systems(OnEnter(State::Menu), setup)
            .add_systems(OnExit(State::Menu), cleanup)
            .add_systems(Update, (run, exit_app).run_if(in_state(State::Menu)));
    }
}

fn setup(mut commands: Commands) {
    let play_button = commands
        .spawn(NodeBundle {
            style: Style {
                // center button
                width: Val::Percent(100.),
                height: Val::Percent(100.),
                justify_content: JustifyContent::Center,
                align_items: AlignItems::Center,
                padding: UiRect::px(10., 10., 10., 10.),
                ..default()
            },
            ..default()
        })
        .with_children(|parent| {
            parent
                .spawn((
                    ButtonBundle {
                        style: Style {
                            width: Val::Px(150.),
                            height: Val::Px(65.),
                            // horizontally center child text
                            justify_content: JustifyContent::Center,
                            // vertically center child text
                            align_items: AlignItems::Center,
                            ..default()
                        },
                        background_color: NORMAL_BUTTON.into(),
                        ..default()
                    },
                    ButtonAction {
                        action: ButtonActions::Play,
                    },
                ))
                .with_children(|parent| {
                    parent.spawn(TextBundle::from_section(
                        "Play",
                        TextStyle {
                            font_size: 40.0,
                            color: Color::rgb(0.9, 0.9, 0.9),
                            ..default()
                        },
                    ));
                });
        })
        .id();

    let quit_button = commands
        .spawn(NodeBundle {
            style: Style {
                // center button
                width: Val::Percent(100.),
                height: Val::Percent(100.),
                justify_content: JustifyContent::Center,
                align_items: AlignItems::Center,
                padding: UiRect::px(10., 10., 10., 10.),
                ..default()
            },
            ..default()
        })
        .with_children(|parent| {
            parent
                .spawn((
                    ButtonBundle {
                        style: Style {
                            width: Val::Px(150.),
                            height: Val::Px(65.),
                            // horizontally center child text
                            justify_content: JustifyContent::Center,
                            // vertically center child text
                            align_items: AlignItems::Center,
                            ..default()
                        },
                        background_color: NORMAL_BUTTON.into(),
                        ..default()
                    },
                    ButtonAction {
                        action: ButtonActions::Quit,
                    },
                ))
                .with_children(|parent| {
                    parent.spawn(TextBundle::from_section(
                        "Quit",
                        TextStyle {
                            font_size: 40.0,
                            color: Color::rgb(0.9, 0.9, 0.9),
                            ..default()
                        },
                    ));
                });
        })
        .id();

    let mut window = commands.spawn(NodeBundle {
        style: Style {
            // center panel
            justify_content: JustifyContent::Center,
            flex_direction: FlexDirection::Column,
            flex_grow: 1.0,
            align_items: AlignItems::Center,
            align_self: AlignSelf::Center,
            justify_self: JustifySelf::Center,

            ..default()
        },
        ..default()
    });

    window.add_child(play_button);
    window.add_child(quit_button);
    let root = window.id();

    commands.insert_resource(MenuData { root });
}

fn run(
    mut next_state: ResMut<NextState<State>>,
    mut app_exit: EventWriter<AppExit>,
    mut interaction_query: Query<
        (&Interaction, &mut BackgroundColor, &ButtonAction),
        (Changed<Interaction>, With<Button>),
    >,
) {
    for (interaction, mut color, action) in &mut interaction_query {
        match *interaction {
            Interaction::Pressed => {
                *color = PRESSED_BUTTON.into();

                match action.action {
                    ButtonActions::Play => next_state.set(State::Game),
                    ButtonActions::Quit => _=app_exit.send(AppExit),
                }
            }
            Interaction::Hovered => {
                *color = HOVERED_BUTTON.into();
            }
            Interaction::None => {
                *color = NORMAL_BUTTON.into();
            }
        }
    }
}

fn cleanup(mut commands: Commands, menu_data: Res<MenuData>) {
    commands.entity(menu_data.root).despawn_recursive();
}

fn exit_app(
    keyboard: Res<ButtonInput<KeyCode>>,
    mut app_exit: EventWriter<AppExit>,
    mut next_state: ResMut<NextState<State>>,
) {
    if keyboard.just_pressed(KeyCode::Escape) {
        app_exit.send(AppExit);
    }

    if keyboard.just_pressed(KeyCode::Space) || keyboard.just_pressed(KeyCode::Enter) {
        next_state.set(State::Game);
    }
}
