/*
 * SPDX-FileCopyrightText: 2026 CSC – IT Center for Science
 *
 * SPDX-License-Identifier: Apache-2.0
 */

#pragma once

#define HPC_CONCAT2_(a, b) a##b

/** Pastes two tokens after macro-expanding them. */
#define HPC_CONCAT2(a, b) HPC_CONCAT2_(a, b)

#define HPC_CONCAT3_(a, b, c) a##b##c

/** Pastes three tokens after macro-expanding them. */
#define HPC_CONCAT3(a, b, c) HPC_CONCAT3_(a, b, c)
