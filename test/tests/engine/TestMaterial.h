/********************************************************************************
* ReactPhysics3D physics library, http://www.reactphysics3d.com                 *
* Copyright (c) 2010-2024 Daniel Chappuis                                       *
*********************************************************************************
*                                                                               *
* This software is provided 'as-is', without any express or implied warranty.   *
* In no event will the authors be held liable for any damages arising from the  *
* use of this software.                                                         *
*                                                                               *
* Permission is granted to anyone to use this software for any purpose,         *
* including commercial applications, and to alter it and redistribute it        *
* freely, subject to the following restrictions:                                *
*                                                                               *
* 1. The origin of this software must not be misrepresented; you must not claim *
*    that you wrote the original software. If you use this software in a        *
*    product, an acknowledgment in the product documentation would be           *
*    appreciated but is not required.                                           *
*                                                                               *
* 2. Altered source versions must be plainly marked as such, and must not be    *
*    misrepresented as the original software.                                     *
* 3. This notice may not be removed or altered from any source distribution.    *
*                                                                               *
********************************************************************************/

#ifndef TEST_MATERIAL_H
#define TEST_MATERIAL_H

// Libraries
#include <reactphysics3d/reactphysics3d.h>
#include <vector>

/// Reactphysics3D namespace
namespace reactphysics3d {

// Class TestMaterial
/**
 * Unit tests for runtime material changes. A material can be changed while the
 * simulation is running (for instance when the ground material is replaced).
 * The new friction coefficient and bounciness must be taken into account at the
 * next solving stage and all the consumers (contact solver, contact caches and
 * queries) must read the same material version. In particular, the warm-starting
 * impulses cached in the contact manifolds must not keep using an outdated
 * material snapshot for more than one committed step.
 */
class TestMaterial : public Test {

    private :

        // ---------- Attributes ---------- //

        PhysicsCommon mPhysicsCommon;

        /// Fixed time step used by the tests
        const decimal mTimeStep;

    public :

        // ---------- Methods ---------- //

        /// Constructor
        TestMaterial(const std::string& name)
            : Test(name), mTimeStep(decimal(1.0 / 60.0)) {

        }

        /// Destructor
        virtual ~TestMaterial() override = default;

        /// Run the tests
        void run() override {
            testMaterialVersion();
            testRuntimeFrictionChange();
            testRuntimeBouncinessChange();
            testStaleContactCacheInvalidation();
            testSleepingBodyUsesNewMaterial();
        }

        /// Create a static ground body with a single box collider
        RigidBody* createGround(PhysicsWorld* world, Collider*& collider,
                                decimal friction, decimal bounciness) {

            RigidBody* ground = world->createRigidBody(Transform(Vector3(0, decimal(-0.5), 0), Quaternion::identity()));
            ground->setType(BodyType::STATIC);
            BoxShape* groundShape = mPhysicsCommon.createBoxShape(Vector3(30, decimal(0.5), 30));
            collider = ground->addCollider(groundShape, Transform::identity());
            collider->getMaterial().setFrictionCoefficient(friction);
            collider->getMaterial().setBounciness(bounciness);
            return ground;
        }

        /// Create a dynamic unit box (mass 1 for the default mass density)
        RigidBody* createBox(PhysicsWorld* world, const Vector3& position, Collider*& collider,
                             decimal friction, decimal bounciness) {

            RigidBody* body = world->createRigidBody(Transform(position, Quaternion::identity()));
            BoxShape* boxShape = mPhysicsCommon.createBoxShape(Vector3(decimal(0.5), decimal(0.5), decimal(0.5)));
            collider = body->addCollider(boxShape, Transform::identity());
            collider->getMaterial().setFrictionCoefficient(friction);
            collider->getMaterial().setBounciness(bounciness);
            return body;
        }

        /// The material version is a single source of truth that is incremented
        /// each time a contact parameter (friction or bounciness) is modified.
        void testMaterialVersion() {

            PhysicsWorld* world = mPhysicsCommon.createPhysicsWorld();
            Collider* groundCollider;
            createGround(world, groundCollider, decimal(0.3), decimal(0.5));

            const uint64 initialVersion = groundCollider->getMaterial().getMaterialVersion();
            rp3d_test(initialVersion > 0);

            // Changing the friction publishes a new version
            groundCollider->getMaterial().setFrictionCoefficient(decimal(0.42));
            rp3d_test(groundCollider->getMaterial().getMaterialVersion() == initialVersion + 1);
            rp3d_test(approxEqual(groundCollider->getMaterial().getFrictionCoefficient(), decimal(0.42)));

            // Changing the bounciness also publishes a new version
            groundCollider->getMaterial().setBounciness(decimal(0.31));
            rp3d_test(groundCollider->getMaterial().getMaterialVersion() == initialVersion + 2);
            rp3d_test(approxEqual(groundCollider->getMaterial().getBounciness(), decimal(0.31)));

            // Setting the same value still publishes a new version (parameter write)
            groundCollider->getMaterial().setBounciness(decimal(0.31));
            rp3d_test(groundCollider->getMaterial().getMaterialVersion() == initialVersion + 3);

            // The query path always reads the latest published values
            rp3d_test(approxEqual(groundCollider->getMaterial().getFrictionCoefficient(), decimal(0.42)));
            rp3d_test(approxEqual(groundCollider->getMaterial().getBounciness(), decimal(0.31)));

            mPhysicsCommon.destroyPhysicsWorld(world);
        }

        /// A runtime friction change must be taken into account at the next solving stage.
        /// The body keeps sliding at constant speed once the ground becomes frictionless.
        void testRuntimeFrictionChange() {

            PhysicsWorld* world = mPhysicsCommon.createPhysicsWorld();
            Collider* groundCollider;
            createGround(world, groundCollider, decimal(0.8), decimal(0.0));

            Collider* boxCollider;
            RigidBody* box = createBox(world, Vector3(0, decimal(0.5), 0), boxCollider, decimal(0.8), decimal(0.0));
            box->setLinearVelocity(Vector3(3, 0, 0));

            // Settle the contact vertically and build up the warm-starting caches
            for (int i = 0; i < 10; i++) {
                world->update(mTimeStep);
            }

            const decimal speedBeforeChange = box->getLinearVelocity().x;
            rp3d_test(speedBeforeChange < decimal(3.0));

            // Replace the ground material with a frictionless one at runtime
            groundCollider->getMaterial().setFrictionCoefficient(decimal(0.0));

            // From the next committed step, the speed must not decrease anymore
            world->update(mTimeStep);
            const decimal speedAtSwitchStep = box->getLinearVelocity().x;
            for (int i = 0; i < 10; i++) {
                world->update(mTimeStep);
            }
            const decimal speedAfterSwitch = box->getLinearVelocity().x;

            rp3d_test(approxEqual(speedAtSwitchStep, speedBeforeChange, decimal(1e-6)));
            rp3d_test(approxEqual(speedAfterSwitch, speedAtSwitchStep, decimal(1e-6)));

            // The query path reads the same material value as the solver
            rp3d_test(approxEqual(groundCollider->getMaterial().getFrictionCoefficient(), decimal(0.0)));

            mPhysicsCommon.destroyPhysicsWorld(world);
        }

        /// A runtime bounciness change must be used for the next impacts.
        void testRuntimeBouncinessChange() {

            PhysicsWorld* world = mPhysicsCommon.createPhysicsWorld();
            Collider* groundCollider;
            createGround(world, groundCollider, decimal(0.0), decimal(0.0));

            // The box starts above the ground and falls on it after the material change
            Collider* boxCollider;
            RigidBody* box = createBox(world, Vector3(0, 3, 0), boxCollider, decimal(0.0), decimal(0.0));

            // Make the ground bouncy at runtime before the first impact
            groundCollider->getMaterial().setBounciness(decimal(0.6));

            // Measure the impact and the rebound velocities
            decimal minNormalVelocity = decimal(0.0);
            decimal maxReboundVelocity = decimal(0.0);
            for (int i = 0; i < 90; i++) {
                world->update(mTimeStep);
                const decimal normalVelocity = box->getLinearVelocity().y;
                if (normalVelocity < minNormalVelocity) minNormalVelocity = normalVelocity;
                if (normalVelocity > maxReboundVelocity) maxReboundVelocity = normalVelocity;
            }

            // The impact is significant and the rebound ratio matches the new bounciness (0.6)
            rp3d_test(minNormalVelocity < decimal(-2.0));
            rp3d_test(maxReboundVelocity > decimal(1.0));
            const decimal reboundRatio = maxReboundVelocity / (-minNormalVelocity);
            rp3d_test(reboundRatio > decimal(0.55));
            rp3d_test(reboundRatio < decimal(0.65));

            mPhysicsCommon.destroyPhysicsWorld(world);
        }

        /// Warm-starting impulses cached in a persistent contact manifold are derived
        /// from the material in use when they are computed. After a material change,
        /// the first committed step must solve the contact as if it were a new contact
        /// (the cached impulses derived from the old material are discarded).
        void testStaleContactCacheInvalidation() {

            // Reference world: the contact is directly created with the target friction
            // coefficient, so no warm-starting impulse computed with another friction
            // coefficient exists.
            PhysicsCommon referenceCommon;
            PhysicsWorld* referenceWorld = referenceCommon.createPhysicsWorld();

            auto velocityChangeAtTargetFriction = [&](PhysicsCommon& common, PhysicsWorld* physicsWorld,
                                                      decimal initialGroundFriction,
                                                      decimal targetGroundFriction) {

                Collider* groundCollider;
                RigidBody* ground = physicsWorld->createRigidBody(
                    Transform(Vector3(0, decimal(-0.5), 0), Quaternion::identity()));
                ground->setType(BodyType::STATIC);
                BoxShape* groundShape = common.createBoxShape(Vector3(30, decimal(0.5), 30));
                groundCollider = ground->addCollider(groundShape, Transform::identity());
                groundCollider->getMaterial().setFrictionCoefficient(initialGroundFriction);

                RigidBody* body = physicsWorld->createRigidBody(
                    Transform(Vector3(0, decimal(0.5), 0), Quaternion::identity()));
                BoxShape* boxShape = common.createBoxShape(
                    Vector3(decimal(0.5), decimal(0.5), decimal(0.5)));
                Collider* boxCollider = body->addCollider(boxShape, Transform::identity());
                boxCollider->getMaterial().setFrictionCoefficient(decimal(0.8));
                body->setLinearVelocity(Vector3(1, 0, 0));

                // Build the persistent contact and its warm-starting caches
                for (int i = 0; i < 5; i++) {
                    physicsWorld->update(mTimeStep);
                }

                const decimal velocityBeforeTargetStep = body->getLinearVelocity().x;
                if (initialGroundFriction != targetGroundFriction) {
                    groundCollider->getMaterial().setFrictionCoefficient(targetGroundFriction);
                }

                // Committed step solved with the target friction coefficient
                physicsWorld->update(mTimeStep);
                return body->getLinearVelocity().x - velocityBeforeTargetStep;
            };

            // Cold start: contact created directly with friction 0.4
            const decimal coldStartVelocityChange =
                velocityChangeAtTargetFriction(referenceCommon, referenceWorld,
                                               decimal(0.4), decimal(0.4));

            // Warm cache: contact history built with friction 0.8, then switched to 0.4.
            // With the stale cache invalidation, the switch step must match the cold start.
            PhysicsWorld* world = mPhysicsCommon.createPhysicsWorld();
            const decimal warmCacheVelocityChange =
                velocityChangeAtTargetFriction(mPhysicsCommon, world,
                                               decimal(0.8), decimal(0.4));

            rp3d_test(approxEqual(warmCacheVelocityChange, coldStartVelocityChange, decimal(2e-4)));

            mPhysicsCommon.destroyPhysicsWorld(world);
            referenceCommon.destroyPhysicsWorld(referenceWorld);
        }

        /// A sleeping body does not need any material update while asleep but the next
        /// contact that wakes it up must be solved with the current material version.
        void testSleepingBodyUsesNewMaterial() {

            PhysicsWorld* world = mPhysicsCommon.createPhysicsWorld();
            Collider* groundCollider;
            createGround(world, groundCollider, decimal(0.8), decimal(0.0));

            Collider* sleeperCollider;
            RigidBody* sleeper = createBox(world, Vector3(0, decimal(0.5), 0),
                                           sleeperCollider, decimal(0.8), decimal(0.0));

            // Wait long enough for the body to fall asleep (default time before sleep is 1s)
            for (int i = 0; i < 120; i++) {
                world->update(mTimeStep);
            }
            rp3d_test(sleeper->isSleeping());

            // Change the ground material while the body is sleeping
            groundCollider->getMaterial().setFrictionCoefficient(decimal(0.0));
            groundCollider->getMaterial().setBounciness(decimal(0.6));
            rp3d_test(sleeper->isSleeping());

            // The next contact (incoming dynamic body) wakes the sleeping body up
            Collider* bulletCollider;
            RigidBody* bullet = createBox(world, Vector3(0, decimal(0.5), decimal(-5.0)),
                                          bulletCollider, decimal(0.0), decimal(0.0));
            bullet->setLinearVelocity(Vector3(0, 0, 8));

            bool hasBeenWokenUp = false;
            for (int i = 0; i < 60; i++) {
                world->update(mTimeStep);
                if (!sleeper->isSleeping()) {
                    hasBeenWokenUp = true;
                    break;
                }
            }
            rp3d_test(hasBeenWokenUp);

            // Once awake, a fresh impact against the ground uses the new bounciness (0.6).
            // Move the body above the ground to get an impact uncontaminated by the
            // split-impulse position correction of the previous resting contact.
            sleeper->setIsSleeping(false);
            sleeper->setTransform(Transform(Vector3(5, 3, 5), Quaternion::identity()));
            sleeper->setLinearVelocity(Vector3(0, 0, 0));
            for (int i = 0; i < 5; i++) {
                world->update(mTimeStep);
            }

            decimal minNormalVelocity = decimal(0.0);
            decimal maxReboundVelocity = decimal(0.0);
            for (int i = 0; i < 90; i++) {
                world->update(mTimeStep);
                const decimal normalVelocity = sleeper->getLinearVelocity().y;
                if (normalVelocity < minNormalVelocity) minNormalVelocity = normalVelocity;
                if (normalVelocity > maxReboundVelocity) maxReboundVelocity = normalVelocity;
            }

            rp3d_test(minNormalVelocity < decimal(-2.0));
            const decimal reboundRatio = maxReboundVelocity / (-minNormalVelocity);
            rp3d_test(reboundRatio > decimal(0.55));
            rp3d_test(reboundRatio < decimal(0.65));

            mPhysicsCommon.destroyPhysicsWorld(world);
        }
};

}

#endif
