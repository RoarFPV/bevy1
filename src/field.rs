use bevy::prelude::*;

#[derive(Resource, Debug, Default)]
pub struct FieldConfig {
    pub size: UVec3,
    pub grid_size: f32,
}

impl FieldConfig {
    fn default() -> Self {
        Self {
            size: UVec3::new(10, 10, 1),
            grid_size: 1.0,
        }
    }
}

#[derive(Resource, Debug)]
pub struct Field {
    pub cells: Vec<Entity>,
    pub size: IVec3,
}

impl Field {
    fn cell_index(&self, point: IVec3) -> i32 {
        let safep = point.clamp(IVec3::ZERO, self.size.saturating_sub(IVec3::ONE));
        safep.x + safep.y * self.size.x as i32 + safep.z * self.size.x as i32 * self.size.y as i32
    }

    pub fn get_entity(&self, point: IVec3) -> Entity {
        let index = self.cell_index(point) as usize;
        self.cells[index]
    }

    pub fn set_entity(&mut self, point: IVec3, entity: Entity) {
        let index = self.cell_index(point) as usize;
        self.cells[index] = entity;
    }

    pub fn new(size: IVec3) -> Field {
        Field {
            cells: vec![Entity::PLACEHOLDER; (size.x * size.y * size.z) as usize],
            size,
        }
    }

    pub fn center(&self) -> Vec3 {
        self.size.as_vec3() / 2.0
    }
}
