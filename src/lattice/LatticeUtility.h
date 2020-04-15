

namespace OML {

	// This does seem to be working..
	struct GLMVec3KeyFuncs
	{
		size_t operator()(const glm::vec3& k) const
		{
			glm::vec3 v;
			v.x = std::roundf(k.x * 1000) / 1000.0f;
			v.y = std::roundf(k.y * 1000) / 1000.0f;
			v.z = std::roundf(k.z * 1000) / 1000.0f;
			return std::hash<glm::vec3>()(v);
		}

		bool operator()(const glm::vec3& a, const glm::vec3& b) const
		{
			return std::abs(a.x - b.x) < 1e-2 &&
				std::abs(a.y - b.y) < 1e-2 &&
				std::abs(a.z - b.z) < 1e-2;
		}
	};
	typedef std::unordered_map<glm::vec3, size_t, GLMVec3KeyFuncs, GLMVec3KeyFuncs> vec3_map;

	/*
		Struct that holds the values to clamp u and v to when evaluating a local surface
		in the shader. In the corners the values would be us = 0.0, ue = 1.0, vs = 0.0 and ve = 1.0

		TODO: Change so patches can reuse BoundaryInfos, they already contain an index to they
		should be able to point to the same. For a regular rectangular grid most BoudnaryInfos will
		be the same, so a lot of memory is wasted. Check if changing this will change performance also.
	*/
	struct BoundaryInfo
	{
		BoundaryInfo()
			: us(0.0f), ue(1.0f), vs(0.0f), ve(1.0f) {}
		BoundaryInfo(float us, float ue, float vs, float ve)
			: us(us), ue(ue), vs(vs), ve(ve) {}

		float us, ue, vs, ve;

		bool operator==(const BoundaryInfo& other) const {
			return std::abs(this->us - other.us) < 1e-2 &&
				std::abs(this->ue - other.ue) < 1e-2 &&
				std::abs(this->vs - other.vs) < 1e-2 &&
				std::abs(this->ve - other.ve) < 1e-2;
		}
	};

	struct BoundaryKeyFunc
	{
		size_t operator()(const BoundaryInfo& b) const
		{
			glm::vec4 v;
			v.x = std::roundf(b.us * 1000) / 1000.0f;
			v.y = std::roundf(b.ue * 1000) / 1000.0f;
			v.z = std::roundf(b.vs * 1000) / 1000.0f;
			v.w = std::roundf(b.ve * 1000) / 1000.0f;
			return std::hash<glm::vec4>()(v);
		}

		bool operator()(const BoundaryInfo& a, const BoundaryInfo& b) const
		{
			return std::abs(a.us - b.us) < 1e-2 &&
				std::abs(a.ue - b.ue) < 1e-2 &&
				std::abs(a.vs - b.vs) < 1e-2 &&
				std::abs(a.ve - b.ve) < 1e-2;
		}
	};

	typedef std::unordered_map<BoundaryInfo, size_t, BoundaryKeyFunc, BoundaryKeyFunc> boundary_map;

}