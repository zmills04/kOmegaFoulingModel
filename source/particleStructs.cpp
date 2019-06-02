#include "particleStructs.h"
#include "clVariablesTR.h"
//////////////////////////////////////////////////////
///////////             PAR              /////////////
//////////////////////////////////////////////////////


void Par::setStruct(legacyPar& struct_, const int i)
{
	pos(i) = struct_.pos;
	Num_rep(i) = struct_.Num_rep;
	type(i) = struct_.type;
	Dep_Flag(i) = struct_.Dep_Flag;
	Dep_timer(i) = struct_.Dep_timer;
	timer(i) = struct_.timer;
	loc(i) = struct_.loc;
}

legacyPar Par::getStruct(const int i)
{
	legacyPar struct_;
	struct_.pos = pos(i);
	struct_.Num_rep = Num_rep(i);
	struct_.type = type(i);
	struct_.Dep_Flag = Dep_Flag(i);
	struct_.Dep_timer = Dep_timer(i);
	struct_.timer = timer(i);
	struct_.loc = loc(i);
	return struct_;
}

void Par::allocateArrays()
{
	pos.setName("Par_Pos");
	Num_rep.setName("Par_NumRep");
	type.setName("Par_Type");
	Dep_Flag.setName("Par_DepFlag");
	Dep_timer.setName("Par_DepTimer");
	timer.setName("Par_Timer");
	loc.setName("Par_Loc");

	pos.allocate(fullSize);
	pos.fill({ { 0., 0. } });
	Num_rep.zeros(fullSize);
	type.zeros(fullSize);
	Dep_Flag.zeros(fullSize);
	Dep_timer.zeros(fullSize);
	timer.zeros(fullSize);
	loc.zeros(fullSize);
}

void Par::allocateBuffers(int bufSize = -1)
{
	if (bufSize == -1)
		bufSize = fullSize;

	bufFullSize = bufSize;

	pos.allocate_buffer_size(bufFullSize);
	Num_rep.allocate_buffer_size(bufFullSize);
	type.allocate_buffer_size(bufFullSize);
	Dep_Flag.allocate_buffer_size(bufFullSize);
	Dep_timer.allocate_buffer_size(bufFullSize);
	timer.allocate_buffer_size(bufFullSize);
	loc.allocate_buffer_size(bufFullSize);
}

void Par::copyToDevice(int writeSize = -1, cl_command_queue * que_ = nullptr,
	cl_bool bFlag_ = CL_TRUE)
{
	if (writeSize == -1)
		writeSize = fullSize;
	pos.copy_to_buffer_size(writeSize, que_, bFlag_);
	Num_rep.copy_to_buffer_size(writeSize, que_, bFlag_);
	type.copy_to_buffer_size(writeSize, que_, bFlag_);
	Dep_Flag.copy_to_buffer_size(writeSize, que_, bFlag_);
	Dep_timer.copy_to_buffer_size(writeSize, que_, bFlag_);
	timer.copy_to_buffer_size(writeSize, que_, bFlag_);
	loc.copy_to_buffer_size(writeSize, que_, bFlag_);
}

void Par::copyToHost(int readSize = -1, cl_command_queue * que_ = nullptr,
	cl_bool bFlag_ = CL_TRUE)
{
	if (readSize == -1)
		readSize = fullSize;
	pos.read_from_buffer_size(readSize, que_, bFlag_);
	Num_rep.read_from_buffer_size(readSize, que_, bFlag_);
	type.read_from_buffer_size(readSize, que_, bFlag_);
	Dep_Flag.read_from_buffer_size(readSize, que_, bFlag_);
	Dep_timer.read_from_buffer_size(readSize, que_, bFlag_);
	timer.read_from_buffer_size(readSize, que_, bFlag_);
	loc.read_from_buffer_size(readSize, que_, bFlag_);
}

void Par::writeParToBuffer(Par& Ptemp, int writeSize = -1, cl_command_queue* que_ = nullptr,
	cl_bool bFlag_ = CL_TRUE)
{
	pos.write_array_to_buffer(Ptemp.pos.get_array(), bFlag_, writeSize, que_);
	Num_rep.write_array_to_buffer(Ptemp.Num_rep.get_array(), bFlag_, writeSize, que_);
	type.write_array_to_buffer(Ptemp.type.get_array(), bFlag_, writeSize, que_);
	Dep_Flag.write_array_to_buffer(Ptemp.Dep_Flag.get_array(), bFlag_, writeSize, que_);
	Dep_timer.write_array_to_buffer(Ptemp.Dep_timer.get_array(), bFlag_, writeSize, que_);
	timer.write_array_to_buffer(Ptemp.timer.get_array(), bFlag_, writeSize, que_);
	loc.write_array_to_buffer(Ptemp.loc.get_array(), bFlag_, writeSize, que_);
}

void Par::copyToParOnDevice(Par& Ptemp, bool copyLoc, int writeSize,
	cl_command_queue* que_, int num_wait, cl_event* wait, cl_event* evt)
{
	//only first event needs to wait, since all on same queue
	copyToArrayOnDevice(Ptemp, posArr, writeSize, que_, num_wait, wait);
	copyToArrayOnDevice(Ptemp, numRepArr, writeSize, que_);
	copyToArrayOnDevice(Ptemp, typeArr, writeSize, que_);
	copyToArrayOnDevice(Ptemp, depFlagArr, writeSize, que_);
	copyToArrayOnDevice(Ptemp, depTimerArr, writeSize, que_);
	copyToArrayOnDevice(Ptemp, timerArr, writeSize, que_);
	//only final call needs to be tracked by evt since all on same queue
	if (!copyLoc)
	{
		copyToArrayOnDevice(Ptemp, timerArr, writeSize, que_, 0, nullptr, evt);
	}
	else
	{
		copyToArrayOnDevice(Ptemp, timerArr, writeSize, que_);
		copyToArrayOnDevice(Ptemp, locArr, writeSize, que_, 0, nullptr, evt);
	}
}

void Par::copyToArrayOnDevice(Par& Ptemp, arrName arrname_, int writeSize,
	cl_command_queue* que_, int num_wait, cl_event* wait, cl_event* evt)
{
	switch (arrname_)
	{
	case posArr:
	{
		pos.enqueue_copy_to_buffer(Ptemp.pos.get_buffer(),
			writeSize, que_, num_wait, wait, evt);
		break;
	}
	case numRepArr:
	{
		Num_rep.enqueue_copy_to_buffer(Ptemp.Num_rep.get_buffer(),
			writeSize, que_, num_wait, wait, evt);
		break;
	}
	case typeArr:
	{
		type.enqueue_copy_to_buffer(Ptemp.type.get_buffer(),
			writeSize, que_, num_wait, wait, evt);
		break;
	}
	case depFlagArr:
	{
		Dep_Flag.enqueue_copy_to_buffer(Ptemp.Dep_Flag.get_buffer(),
			writeSize, que_, num_wait, wait, evt);
		break;
	}
	case depTimerArr:
	{
		Dep_timer.enqueue_copy_to_buffer(Ptemp.Dep_timer.get_buffer(),
			writeSize, que_, num_wait, wait, evt);
		break;
	}
	case timerArr:
	{
		timer.enqueue_copy_to_buffer(Ptemp.timer.get_buffer(),
			writeSize, que_, num_wait, wait, evt);
		break;
	}
	case locArr:
	{
		loc.enqueue_copy_to_buffer(Ptemp.loc.get_buffer(),
			writeSize, que_, num_wait, wait, evt);
		break;
	}
	}
}

void Par::copyToParOnDeviceBlocking(Par& Ptemp, bool copyLoc, int writeSize,
	cl_command_queue* que_, int num_wait, cl_event* wait)
{
	//only first event needs to wait, since all on same queue
	copyToArrayOnDevice(Ptemp, posArr, writeSize, que_, num_wait, wait);
	copyToArrayOnDevice(Ptemp, numRepArr, writeSize, que_);
	copyToArrayOnDevice(Ptemp, typeArr, writeSize, que_);
	copyToArrayOnDevice(Ptemp, depFlagArr, writeSize, que_);
	copyToArrayOnDevice(Ptemp, depTimerArr, writeSize, que_);

	//only last call needs to be blocking because all in same queue
	if (!copyLoc)
	{
		copyToArrayOnDeviceBlocking(Ptemp, timerArr, writeSize, que_);
	}
	else
	{
		copyToArrayOnDevice(Ptemp, timerArr, writeSize, que_);
		copyToArrayOnDeviceBlocking(Ptemp, locArr, writeSize, que_);
	}
}

void Par::copyToArrayOnDeviceBlocking(Par& Ptemp, arrName arrname_, int writeSize,
	cl_command_queue* que_, int num_wait, cl_event* wait)
{
	switch (arrname_)
	{
	case posArr:
	{
		pos.enqueue_copy_to_buffer_blocking(Ptemp.pos.get_buffer(),
			writeSize, que_, num_wait, wait);
		break;
	}
	case numRepArr:
	{
		Num_rep.enqueue_copy_to_buffer_blocking(Ptemp.Num_rep.get_buffer(),
			writeSize, que_, num_wait, wait);
		break;
	}
	case typeArr:
	{
		type.enqueue_copy_to_buffer_blocking(Ptemp.type.get_buffer(),
			writeSize, que_, num_wait, wait);
		break;
	}
	case depFlagArr:
	{
		Dep_Flag.enqueue_copy_to_buffer_blocking(Ptemp.Dep_Flag.get_buffer(),
			writeSize, que_, num_wait, wait);
		break;
	}
	case depTimerArr:
	{
		Dep_timer.enqueue_copy_to_buffer_blocking(Ptemp.Dep_timer.get_buffer(),
			writeSize, que_, num_wait, wait);
		break;
	}
	case timerArr:
	{
		timer.enqueue_copy_to_buffer_blocking(Ptemp.timer.get_buffer(),
			writeSize, que_, num_wait, wait);
		break;
	}
	case locArr:
	{
		loc.enqueue_copy_to_buffer_blocking(Ptemp.loc.get_buffer(),
			writeSize, que_, num_wait, wait);
		break;
	}
	}
}


void Par::setBuffers(Kernel& ker, int& curind, arrName arrList[], int numArrs)
{
	for (int i = 0; i < numArrs; i++)
	{
		switch (arrList[i])
		{
		case posArr:
		{
			ker.set_argument(curind, pos.get_buf_add());
			break;
		}
		case numRepArr:
		{
			ker.set_argument(curind, Num_rep.get_buf_add());
			break;
		}
		case typeArr:
		{
			ker.set_argument(curind, type.get_buf_add());
			break;
		}
		case depFlagArr:
		{
			ker.set_argument(curind, Dep_Flag.get_buf_add());
			break;
		}
		case depTimerArr:
		{
			ker.set_argument(curind, Dep_timer.get_buf_add());
			break;
		}
		case timerArr:
		{
			ker.set_argument(curind, timer.get_buf_add());
			break;
		}
		case locArr:
		{
			ker.set_argument(curind, loc.get_buf_add());
			break;
		}
		}
		curind++;
	}
}

void Par::ini()
{
	//initialization functions specific to particles
}

bool Par::load()
{
	int ret = true;
	ret &= pos.load("load" SLASH + pos.getName());
	ret &= Num_rep.load("load" SLASH + Num_rep.getName());
	ret &= type.load("load" SLASH + type.getName());
	ret &= Dep_Flag.load("load" SLASH + Dep_Flag.getName());
	ret &= Dep_timer.load("load" SLASH + Dep_timer.getName());
	ret &= timer.load("load" SLASH + timer.getName());
	return ret;
	for (int i = 0; i < fullSize; i++)
	{
		cl_int2 Posi = { { (int)floor(pos(i).x), (int)floor(pos(i).y) } };
		loc(i) = Posi.x + vtr.TrDomainSize.x*Posi.y;
	}
}

bool Par::save(bool saveAll = false, saveFlags saveFlag_ = saveFl)
{
	bool ret = true;
	switch (saveFlag_)
	{
	case saveFl:
	{
		ret &= pos.save();
		ret &= Num_rep.save();
		ret &= type.save();
		break;
	}
	case saveTxtFl:
	{
		ret &= pos.save2file();
		ret &= Num_rep.save2file();
		ret &= type.save2file();
		break;
	}
	case saveBinFl:
	{
		ret &= pos.savebin();
		ret &= Num_rep.savebin();
		ret &= type.savebin();
		ret &= Dep_Flag.savebin();
		ret &= Dep_timer.savebin();
		ret &= timer.savebin();
		break;
	}
	}

	if (saveAll)
	{
		switch (saveFlag_)
		{
		case saveFl:
		{
			ret &= Dep_Flag.save();
			ret &= Dep_timer.save();
			ret &= timer.save();
			break;
		}
		case saveTxtFl:
		{
			ret &= Dep_Flag.save2file();
			ret &= Dep_timer.save2file();
			ret &= timer.save2file();
			break;
		}
		case saveBinFl:
		{
			break;
		}
		}
	}
	return ret;
}

bool Par::saveFromDevice(bool saveAll = false, saveFlags saveFlag_ = saveFl,
	cl_command_queue * que_ = nullptr)
{
	bool ret = true;
	switch (saveFlag_)
	{
	case saveFl:
	{
		ret &= pos.save_from_device();
		ret &= Num_rep.save_from_device();
		ret &= type.save_from_device();
		break;
	}
	case saveTxtFl:
	{
		ret &= pos.save_txt_from_device("", que_);
		ret &= Num_rep.save_txt_from_device("", que_);
		ret &= type.save_txt_from_device("", que_);
		break;
	}
	case saveBinFl:
	{
		ret &= pos.save_bin_from_device("", que_);
		ret &= Num_rep.save_bin_from_device("", que_);
		ret &= type.save_bin_from_device("", que_);
		ret &= Dep_Flag.save_bin_from_device("", que_);
		ret &= Dep_timer.save_bin_from_device("", que_);
		ret &= timer.save_bin_from_device("", que_);
		break;
	}
	}

	if (saveAll)
	{
		switch (saveFlag_)
		{
		case saveFl:
		{
			ret &= Dep_Flag.save_from_device();
			ret &= Dep_timer.save_from_device();
			ret &= timer.save_from_device();
			break;
		}
		case saveTxtFl:
		{
			ret &= Dep_Flag.save_txt_from_device("", que_);
			ret &= Dep_timer.save_txt_from_device("", que_);
			ret &= timer.save_txt_from_device("", que_);
			break;
		}
		case saveBinFl:
		{
			break;
		}
		}
	}
	return ret;
}

//////////////////////////////////////////////////////
///////////            PParam            /////////////
//////////////////////////////////////////////////////

void PParam::setStruct(legacyPParam& struct_, const int i)
{
	Q_A_prime(i) = struct_.Q_A_prime;
	Q_A(i) = struct_.Q_A;
	tau_crit(i) = struct_.tau_crit;
	Dp(i) = struct_.Dp;
	Mp(i) = struct_.Mp;
	Kth(i) = struct_.Kth;
	D_dist(i) = struct_.D_dist;
	L_coeff(i) = struct_.L_coeff;
	D_coeff(i) = struct_.D_coeff;
}

legacyPParam PParam::getStruct(const int i)
{
	legacyPParam struct_;
	struct_.Q_A_prime = Q_A_prime(i);
	struct_.Q_A = Q_A(i);
	struct_.tau_crit = tau_crit(i);
	struct_.Dp = Dp(i);
	struct_.Mp = Mp(i);
	struct_.Kth = Kth(i);
	struct_.D_dist = D_dist(i);
	struct_.L_coeff = L_coeff(i);
	struct_.D_coeff = D_coeff(i);
	return struct_;
}

void PParam::allocateArrays()
{
	Q_A_prime.setName("PParam_Q_A_prime");
	Q_A.setName("PParam_Q_A");
	tau_crit.setName("PParam_TauCrit");
	Dp.setName("PParam_Dp");
	Mp.setName("PParam_Mp");
	Kth.setName("PParam_Kth");
	D_dist.setName("PParam_D_dist");
	L_coeff.setName("PParam_L_coeff");
	D_coeff.setName("PParam_D_coeff");

	Q_A_prime.allocate(fullSize);
	Q_A.allocate(fullSize);
	tau_crit.allocate(fullSize);
	Dp.zeros(fullSize);
	Mp.zeros(fullSize);
	Kth.zeros(fullSize);
	D_dist.zeros(fullSize);
	L_coeff.zeros(fullSize);
	D_coeff.zeros(fullSize);
}

void PParam::allocateBuffers(int bufSize = -1)
{
	if (bufSize == -1)
		bufSize = fullSize;

	bufFullSize = bufSize;

	Q_A_prime.allocate_buffer_size(bufFullSize);
	Q_A.allocate_buffer_size(bufFullSize);
	tau_crit.allocate_buffer_size(bufFullSize);
	Dp.allocate_buffer_size(bufFullSize);
	Mp.allocate_buffer_size(bufFullSize);
	Kth.allocate_buffer_size(bufFullSize);
	D_dist.allocate_buffer_size(bufFullSize);
	L_coeff.allocate_buffer_size(bufFullSize);
	D_coeff.allocate_buffer_size(bufFullSize);

}

void PParam::copyToDevice(int writeSize = -1, cl_command_queue * que_ = nullptr,
	cl_bool bFlag_ = CL_TRUE)
{
	if (writeSize == -1)
		writeSize = fullSize;

	Q_A_prime.copy_to_buffer_size(writeSize, que_, bFlag_);
	Q_A.copy_to_buffer_size(writeSize, que_, bFlag_);
	tau_crit.copy_to_buffer_size(writeSize, que_, bFlag_);
	Dp.copy_to_buffer_size(writeSize, que_, bFlag_);
	Mp.copy_to_buffer_size(writeSize, que_, bFlag_);
	Kth.copy_to_buffer_size(writeSize, que_, bFlag_);
	D_dist.copy_to_buffer_size(writeSize, que_, bFlag_);
	L_coeff.copy_to_buffer_size(writeSize, que_, bFlag_);
	D_coeff.copy_to_buffer_size(writeSize, que_, bFlag_);
}

void PParam::copyToHost(int readSize = -1, cl_command_queue * que_ = nullptr,
	cl_bool bFlag_ = CL_TRUE)
{
	if (readSize == -1)
		readSize = fullSize;
	Q_A_prime.read_from_buffer_size(readSize, que_, bFlag_);
	Q_A.read_from_buffer_size(readSize, que_, bFlag_);
	tau_crit.read_from_buffer_size(readSize, que_, bFlag_);
	Dp.read_from_buffer_size(readSize, que_, bFlag_);
	Mp.read_from_buffer_size(readSize, que_, bFlag_);
	Kth.read_from_buffer_size(readSize, que_, bFlag_);
	D_dist.read_from_buffer_size(readSize, que_, bFlag_);
	L_coeff.read_from_buffer_size(readSize, que_, bFlag_);
	D_coeff.read_from_buffer_size(readSize, que_, bFlag_);
}

void PParam::ini()
{
	//initialization functions specific to class
}

bool PParam::save(bool saveAll = false, saveFlags saveFlag_ = saveFl)
{
	bool ret = true;

	if (!saveAll)
	{
		std::cout << "Warning: to save PParams, set saveAll to true when calling save\n";
	}

	if (saveAll)
	{
		switch (saveFlag_)
		{
		case saveFl:
		{
			ret &= Q_A_prime.save2file();
			ret &= Q_A.save2file();
			ret &= tau_crit.save2file();
			ret &= Dp.save2file();
			ret &= Mp.save2file();
			ret &= Kth.save2file();
			ret &= D_dist.save2file();
			ret &= L_coeff.save2file();
			ret &= D_coeff.save2file();

			break;
		}
		case saveTxtFl:
		{
			ret &= Q_A_prime.save2file();
			ret &= Q_A.save2file();
			ret &= tau_crit.save2file();
			ret &= Dp.save2file();
			ret &= Mp.save2file();
			ret &= Kth.save2file();
			ret &= D_dist.save2file();
			ret &= L_coeff.save2file();
			ret &= D_coeff.save2file();

			break;
		}
		case saveBinFl:
		{
			std::cout << "Warning: PParam not configured to save to bin files\n";
			break;
		}
		}
	}
	return ret;
}

bool PParam::saveFromDevice(bool saveAll = false, saveFlags saveFlag_ = saveFl,
	cl_command_queue * que_ = nullptr)
{
	std::cout << "Warning: Saving PParams stored on host even though saveFromDevice was called\n";

	return save(saveAll, saveFlag_);
}

void PParam::setBuffers(Kernel& ker, int& curind, arrName arrList[], int numArrs)
{
	for (int i = 0; i < numArrs; i++)
	{
		switch (arrList[i])
		{
		case qaPrimeArr:
		{
			ker.set_argument(curind, Q_A_prime.get_buf_add());
			break;
		}
		case qaArr:
		{
			ker.set_argument(curind, Q_A.get_buf_add());
			break;
		}
		case tauCritArr:
		{
			ker.set_argument(curind, tau_crit.get_buf_add());
			break;
		}
		case dpArr:
		{
			ker.set_argument(curind, Dp.get_buf_add());
			break;
		}
		case mpArr:
		{
			ker.set_argument(curind, Mp.get_buf_add());
			break;
		}
		case kthArr:
		{
			ker.set_argument(curind, Kth.get_buf_add());
			break;
		}
		case distArr:
		{
			ker.set_argument(curind, D_dist.get_buf_add());
			break;
		}
		case dCoeffArr:
		{
			ker.set_argument(curind, D_coeff.get_buf_add());
			break;
		}
		case lCoeffArr:
		{
			ker.set_argument(curind, L_coeff.get_buf_add());
			break;
		}
		}
		curind++;
	}
}



//////////////////////////////////////////////////////
///////////            NodeI             /////////////
//////////////////////////////////////////////////////

void NodeI::setStruct(legacyNodeI& struct_, const int i, const int j)
{
	BLind(i, j) = struct_.BLind;
	wallFlag(i, j) = struct_.wallFlag;
}

legacyNodeI NodeI::getStruct(const int i, const int j)
{
	legacyNodeI struct_;
	struct_.BLind = BLind(i, j);
	struct_.wallFlag = wallFlag(i, j);
	return struct_;
}


void NodeI::allocateArrays()
{
	BLind.setName("NodeI_BLind");
	wallFlag.setName("NodeI_wallFlag");
	BLind.zeros(xSize, xSizeFull, ySize, ySize);
	wallFlag.zeros(xSize, xSizeFull, ySize, ySize);
}

void NodeI::allocateBuffers(int bufSize = -1)
{
	if (bufSize == -1)
		bufSize = fullSize;

	bufFullSize = bufSize;

	BLind.allocate_buffer_size(bufFullSize);
	wallFlag.allocate_buffer_size(bufFullSize);
}

void NodeI::copyToDevice(int writeSize = -1, cl_command_queue * que_ = nullptr,
	cl_bool bFlag_ = CL_TRUE)
{
	if (writeSize == -1)
		writeSize = fullSize;
	BLind.copy_to_buffer_size(writeSize, que_, bFlag_);
	wallFlag.copy_to_buffer_size(writeSize, que_, bFlag_);
}

void NodeI::copyToHost(int readSize = -1, cl_command_queue * que_ = nullptr,
	cl_bool bFlag_ = CL_TRUE)
{
	if (readSize == -1)
		readSize = fullSize;
	BLind.read_from_buffer_size(readSize, que_, bFlag_);
	wallFlag.read_from_buffer_size(readSize, que_, bFlag_);
}

void NodeI::ini()
{
	//initialization functions specific to class
}

bool NodeI::save(bool saveAll = false, saveFlags saveFlag_ = saveFl)
{ //nothing needs to be saved when not debugging, so nothing saved unless saveAll = true
	bool ret = true;
	if (saveAll)
	{
		switch (saveFlag_)
		{
		case saveFl:
		{// no saving bin files, so just call save2file
		}
		case saveTxtFl:
		{
			ret &= BLind.save2file();
			ret &= wallFlag.save2file();
			break;
		}
		case saveBinFl:
		{
			break;
		}
		}
	}
	return ret;
}

bool NodeI::saveFromDevice(bool saveAll = false, saveFlags saveFlag_ = saveFl,
	cl_command_queue * que_ = nullptr)
{
	bool ret = true;
	if (saveAll)
	{
		switch (saveFlag_)
		{
		case saveFl:
		{
		}
		case saveTxtFl:
		{
			ret &= BLind.save_txt_from_device("", que_);
			ret &= wallFlag.save_txt_from_device("", que_);
			break;
		}
		case saveBinFl:
		{
			break;
		}
		}
	}
	return ret;
}


//////////////////////////////////////////////////////
///////////            Neighs            /////////////
//////////////////////////////////////////////////////


void Neighs::setStruct(legacyNeighs& struct_, const int i, const int j)
{
	ii00(i, j) = struct_.ii00;
	ii10(i, j) = struct_.ii10;
	ii01(i, j) = struct_.ii01;
	ii11(i, j) = struct_.ii11;
}

legacyNeighs Neighs::getStruct(const int i, const int j)
{
	legacyNeighs struct_;
	struct_.ii00 = ii00(i, j);
	struct_.ii10 = ii10(i, j);
	struct_.ii01 = ii01(i, j);
	struct_.ii11 = ii11(i, j);
	return struct_;
}


void Neighs::allocateArrays()
{
	ii00.setName("Neighs_ii00");
	ii10.setName("Neighs_ii10");
	ii01.setName("Neighs_ii01");
	ii11.setName("Neighs_ii11");

	ii00.zeros(xSize, xSizeFull, ySize, ySize);
	ii01.zeros(xSize, xSizeFull, ySize, ySize);
	ii10.zeros(xSize, xSizeFull, ySize, ySize);
	ii11.zeros(xSize, xSizeFull, ySize, ySize);
}

void Neighs::allocateBuffers(int bufSize = -1)
{
	if (bufSize == -1)
		bufSize = fullSize;

	bufFullSize = bufSize;

	ii00.allocate_buffer_size(bufFullSize);
	ii10.allocate_buffer_size(bufFullSize);
	ii01.allocate_buffer_size(bufFullSize);
	ii11.allocate_buffer_size(bufFullSize);

}

void Neighs::copyToDevice(int writeSize = -1, cl_command_queue * que_ = nullptr,
	cl_bool bFlag_ = CL_TRUE)
{
	if (writeSize == -1)
		writeSize = fullSize;
	ii00.copy_to_buffer_size(writeSize, que_, bFlag_);
	ii01.copy_to_buffer_size(writeSize, que_, bFlag_);
	ii10.copy_to_buffer_size(writeSize, que_, bFlag_);
	ii11.copy_to_buffer_size(writeSize, que_, bFlag_);
}

void Neighs::copyToHost(int readSize = -1, cl_command_queue * que_ = nullptr,
	cl_bool bFlag_ = CL_TRUE)
{
	if (readSize == -1)
		readSize = fullSize;
	ii00.read_from_buffer_size(readSize, que_, bFlag_);
	ii10.read_from_buffer_size(readSize, que_, bFlag_);
	ii01.read_from_buffer_size(readSize, que_, bFlag_);
	ii11.read_from_buffer_size(readSize, que_, bFlag_);
}

void Neighs::ini()
{
	//initialization functions specific to class
}

bool Neighs::save(bool saveAll = false, saveFlags saveFlag_ = saveFl)
{// no save bin, only saves text  files when saveAll = true
	bool ret = true;

	if (saveAll)
	{
		switch (saveFlag_)
		{
		case saveFl:
		{
		}
		case saveTxtFl:
		{
			ret &= ii00.save2file();
			ret &= ii10.save2file();
			ret &= ii01.save2file();
			ret &= ii11.save2file();
			break;
		}
		case saveBinFl:
		{
			break;
		}
		}
	}
	return ret;
}

bool Neighs::saveFromDevice(bool saveAll = false, saveFlags saveFlag_ = saveFl,
	cl_command_queue * que_ = nullptr)
{
	bool ret = true;

	if (saveAll)
	{
		switch (saveFlag_)
		{
		case saveFl:
		{
		}
		case saveTxtFl:
		{
			ret &= ii00.save_txt_from_device("", que_);
			ret &= ii01.save_txt_from_device("", que_);
			ret &= ii10.save_txt_from_device("", que_);
			ret &= ii11.save_txt_from_device("", que_);
			break;
		}
		case saveBinFl:
		{
			break;
		}
		}
	}
	return ret;
}


//////////////////////////////////////////////////////
///////////            TrParam           /////////////
//////////////////////////////////////////////////////
legacyTrParam TrParam::getStruct()
{
	legacyTrParam struct_;
	struct_.Top_location = fVals(0);
	struct_.Bottom_location = fVals(1);
	struct_.umax_val = fVals(2);
	struct_.bval = fVals(3);
	struct_.offset_y = fVals(4);
	return struct_;
}


void TrParam::allocateArrays()
{
	fVals.setName("TrParam");
	fVals.zeros(fullSize);
}

void TrParam::allocateBuffers(int bufSize = -1)
{
	if (bufSize == -1)
		bufSize = fullSize;

	bufFullSize = bufSize;
	fVals.allocate_buffer_size(bufFullSize);
}


void TrParam::copyToDevice(int writeSize = -1, cl_command_queue * que_ = nullptr,
	cl_bool bFlag_ = CL_TRUE)
{
	if (writeSize == -1)
		writeSize = fullSize;
	fVals.copy_to_buffer_size(writeSize, que_, bFlag_);

}

void TrParam::copyToHost(int readSize = -1, cl_command_queue * que_ = nullptr,
	cl_bool bFlag_ = CL_TRUE)
{
	if (readSize == -1)
		readSize = fullSize;
	fVals.read_from_buffer_size(readSize, que_, bFlag_);

}

void TrParam::ini()
{
	//initialization functions specific to class
}

bool TrParam::save(bool saveAll = false, saveFlags saveFlag_ = saveFl)
{ // only save all saves data
	bool ret = true;


	if (saveAll)
	{
		switch (saveFlag_)
		{
		case saveFl:
		{
		}
		case saveTxtFl:
		{
			ret &= fVals.save2file();
			break;
		}
		case saveBinFl:
		{
			break;
		}
		}
	}
	return ret;
}

bool TrParam::saveFromDevice(bool saveAll = false, saveFlags saveFlag_ = saveFl,
	cl_command_queue * que_ = nullptr)
{
	bool ret = true;

	if (saveAll)
	{
		switch (saveFlag_)
		{
		case saveFl:
		{
		}
		case saveTxtFl:
		{
			ret &= fVals.save_txt_from_device("", que_);

			break;
		}
		case saveBinFl:
		{
			break;
		}
		}
	}
	return ret;
}



//////////////////////////////////////////////////////
///////////            NodeC             /////////////
//////////////////////////////////////////////////////


void NodeC::setStruct(legacyNodeC& struct_, const int i, const int j)
{
	CoeffT00(i, j) = struct_.CoeffT00;
	CoeffT01(i, j) = struct_.CoeffT01;
	CoeffT10(i, j) = struct_.CoeffT10;
	CoeffT11(i, j) = struct_.CoeffT11;
	CoeffU00(i, j) = struct_.CoeffU00;
	CoeffU01(i, j) = struct_.CoeffU01;
	CoeffU10(i, j) = struct_.CoeffU10;
	CoeffU11(i, j) = struct_.CoeffU11;
	neigh(i, j) = struct_.neigh;
}

legacyNodeC NodeC::getStruct(const int i, const int j)
{
	legacyNodeC struct_;
	struct_.CoeffT00 = CoeffT00(i, j);
	struct_.CoeffT01 = CoeffT01(i, j);
	struct_.CoeffT10 = CoeffT10(i, j);
	struct_.CoeffT11 = CoeffT11(i, j);
	struct_.CoeffU00 = CoeffU00(i, j);
	struct_.CoeffU01 = CoeffU01(i, j);
	struct_.CoeffU10 = CoeffU10(i, j);
	struct_.CoeffU11 = CoeffU11(i, j);
	struct_.neigh = neigh(i, j);
	return struct_;
}


void NodeC::allocateArrays()
{
	CoeffT00.setName("NodeC_cT00");
	CoeffT00.allocate(xSize, xSizeFull, ySize, ySize);
	CoeffT00.fill({ {0.,0.,0.,0.} });

	CoeffT01.setName("NodeC_cT01");
	CoeffT01.allocate(xSize, xSizeFull, ySize, ySize);
	CoeffT01.fill({ {0.,0.,0.,0.} });

	CoeffT10.setName("NodeC_cT10");
	CoeffT10.allocate(xSize, xSizeFull, ySize, ySize);
	CoeffT10.fill({ {0.,0.,0.,0.} });

	CoeffT11.setName("NodeC_cT11");
	CoeffT11.allocate(xSize, xSizeFull, ySize, ySize);
	CoeffT11.fill({ {0.,0.,0.,0.} });

	CoeffU00.setName("NodeC_cU00");
	CoeffU00.allocate(xSize, xSizeFull, ySize, ySize);
	CoeffU00.fill({ {0.,0.,0.,0.} });

	CoeffU01.setName("NodeC_cU01");
	CoeffU01.allocate(xSize, xSizeFull, ySize, ySize);
	CoeffU01.fill({ {0.,0.,0.,0.} });

	CoeffU10.setName("NodeC_cU10");
	CoeffU10.allocate(xSize, xSizeFull, ySize, ySize);
	CoeffU10.fill({ {0.,0.,0.,0.} });

	CoeffU11.setName("NodeC_cU11");
	CoeffU11.allocate(xSize, xSizeFull, ySize, ySize);
	CoeffU11.fill({ {0.,0.,0.,0.} });

	neigh.setName("NodeC_neigh");
	neigh.allocate(xSize, xSizeFull, ySize, ySize);
	neigh.fill({ {0,0,0,0} });
}

void NodeC::allocateBuffers(int bufSize = -1)
{
	if (bufSize == -1)
		bufSize = fullSize;

	bufFullSize = bufSize;

	CoeffT00.allocate_buffer_size(bufFullSize);
	CoeffT01.allocate_buffer_size(bufFullSize);
	CoeffT10.allocate_buffer_size(bufFullSize);
	CoeffT11.allocate_buffer_size(bufFullSize);
	CoeffU00.allocate_buffer_size(bufFullSize);
	CoeffU01.allocate_buffer_size(bufFullSize);
	CoeffU10.allocate_buffer_size(bufFullSize);
	CoeffU11.allocate_buffer_size(bufFullSize);
	neigh.allocate_buffer_size(bufFullSize);
}

void NodeC::copyToDevice(int writeSize = -1, cl_command_queue * que_ = nullptr, cl_bool bFlag_ = CL_TRUE)
{
	if (writeSize == -1)
		writeSize = fullSize;

	CoeffT00.copy_to_buffer_size(writeSize, que_, bFlag_);
	CoeffT01.copy_to_buffer_size(writeSize, que_, bFlag_);
	CoeffT10.copy_to_buffer_size(writeSize, que_, bFlag_);
	CoeffT11.copy_to_buffer_size(writeSize, que_, bFlag_);
	CoeffU00.copy_to_buffer_size(writeSize, que_, bFlag_);
	CoeffU01.copy_to_buffer_size(writeSize, que_, bFlag_);
	CoeffU10.copy_to_buffer_size(writeSize, que_, bFlag_);
	CoeffU11.copy_to_buffer_size(writeSize, que_, bFlag_);
	neigh.copy_to_buffer_size(writeSize, que_, bFlag_);
}

void NodeC::copyToHost(int readSize = -1, cl_command_queue * que_ = nullptr, cl_bool bFlag_ = CL_TRUE)
{
	if (readSize == -1)
		readSize = fullSize;


	CoeffT00.read_from_buffer_size(readSize, que_, bFlag_);
	CoeffT01.read_from_buffer_size(readSize, que_, bFlag_);
	CoeffT10.read_from_buffer_size(readSize, que_, bFlag_);
	CoeffT11.read_from_buffer_size(readSize, que_, bFlag_);
	CoeffU00.read_from_buffer_size(readSize, que_, bFlag_);
	CoeffU01.read_from_buffer_size(readSize, que_, bFlag_);
	CoeffU10.read_from_buffer_size(readSize, que_, bFlag_);
	CoeffU11.read_from_buffer_size(readSize, que_, bFlag_);
	neigh.read_from_buffer_size(readSize, que_, bFlag_);
}

void NodeC::ini()
{
	//initialization functions specific to class
}

bool NodeC::save(bool saveAll = false, saveFlags saveFlag_ = saveFl)
{// only need to save for debugging, so only save when saveAll is called
	bool ret = true;
	if (saveAll)
	{
		switch (saveFlag_)
		{
		case saveFl:
		{
		}
		case saveTxtFl:
		{
			ret &= CoeffT00.save2file();
			ret &= CoeffT01.save2file();
			ret &= CoeffT10.save2file();
			ret &= CoeffT11.save2file();
			ret &= CoeffU00.save2file();
			ret &= CoeffU01.save2file();
			ret &= CoeffU10.save2file();
			ret &= CoeffU11.save2file();
			ret &= neigh.save2file();
			break;
		}
		case saveBinFl:
		{
			break;
		}
		}
	}
	return ret;
}

bool NodeC::saveFromDevice(bool saveAll = false, saveFlags saveFlag_ = saveFl,
	cl_command_queue * que_ = nullptr)
{
	bool ret = true;

	if (saveAll)
	{
		switch (saveFlag_)
		{
		case saveFl:
		{
		}
		case saveTxtFl:
		{
			ret &= CoeffT00.save_txt_from_device("", que_);
			ret &= CoeffT01.save_txt_from_device("", que_);
			ret &= CoeffT10.save_txt_from_device("", que_);
			ret &= CoeffT11.save_txt_from_device("", que_);
			ret &= CoeffU00.save_txt_from_device("", que_);
			ret &= CoeffU01.save_txt_from_device("", que_);
			ret &= CoeffU10.save_txt_from_device("", que_);
			ret &= CoeffU11.save_txt_from_device("", que_);
			ret &= neigh.save_txt_from_device("", que_);
			break;
		}
		case saveBinFl:
		{
			break;
		}
		}
	}
	return ret;
}

//////////////////////////////////////////////////////
///////////            NodeV             /////////////
//////////////////////////////////////////////////////
void NodeV::setStruct(legacyNodeV& struct_, const int i, const int j)
{
	Temps(i, j) = struct_.Temps;
	U00(i, j) = struct_.U00;
	U01(i, j) = struct_.U01;
	U10(i, j) = struct_.U10;
	U11(i, j) = struct_.U11;
}

legacyNodeV NodeV::getStruct(const int i, const int j)
{
	legacyNodeV struct_;
	struct_.Temps = Temps(i, j);
	struct_.U00 = U00(i, j);
	struct_.U01 = U01(i, j);
	struct_.U10 = U10(i, j);
	struct_.U11 = U11(i, j);
	return struct_;
}


void NodeV::allocateArrays()
{
	Temps.setName("NodeV_Temps");
	Temps.allocate(xSize, xSizeFull, ySize, ySize);
	Temps.fill({ {0.,0.,0.,0.} });

	U00.setName("NodeV_U00");
	U00.allocate(xSize, xSizeFull, ySize, ySize);
	U00.fill({ {0.,0.} });

	U01.setName("NodeV_U01");
	U01.allocate(xSize, xSizeFull, ySize, ySize);
	U01.fill({ {0.,0.} });

	U10.setName("NodeV_U10");
	U10.allocate(xSize, xSizeFull, ySize, ySize);
	U10.fill({ {0.,0.} });

	U11.setName("NodeV_U11");
	U11.allocate(xSize, xSizeFull, ySize, ySize);
	U11.fill({ {0.,0.} });
}

void NodeV::allocateBuffers(int bufSize = -1)
{
	if (bufSize == -1)
		bufSize = fullSize;

	bufFullSize = bufSize;

	Temps.allocate_buffer_size(bufFullSize);
	U00.allocate_buffer_size(bufFullSize);
	U01.allocate_buffer_size(bufFullSize);
	U10.allocate_buffer_size(bufFullSize);
	U11.allocate_buffer_size(bufFullSize);
}

void NodeV::copyToDevice(int writeSize = -1, cl_command_queue * que_ = nullptr,
	cl_bool bFlag_ = CL_TRUE)
{
	if (writeSize == -1)
		writeSize = fullSize;
	Temps.copy_to_buffer_size(writeSize, que_, bFlag_);
	U00.copy_to_buffer_size(writeSize, que_, bFlag_);
	U01.copy_to_buffer_size(writeSize, que_, bFlag_);
	U10.copy_to_buffer_size(writeSize, que_, bFlag_);
	U11.copy_to_buffer_size(writeSize, que_, bFlag_);
}

void NodeV::copyToHost(int readSize = -1, cl_command_queue * que_ = nullptr, 
	cl_bool bFlag_ = CL_TRUE)
{
	if (readSize == -1)
		readSize = fullSize;
	Temps.read_from_buffer_size(readSize, que_, bFlag_);
	U00.read_from_buffer_size(readSize, que_, bFlag_);
	U01.read_from_buffer_size(readSize, que_, bFlag_);
	U10.read_from_buffer_size(readSize, que_, bFlag_);
	U11.read_from_buffer_size(readSize, que_, bFlag_);
}

void NodeV::ini()
{
	//initialization functions specific to class
}

bool NodeV::save(bool saveAll = false, saveFlags saveFlag_ = saveFl)
{
	bool ret = true;
	if (saveAll)
	{
		switch (saveFlag_)
		{
		case saveFl:
		{
		}
		case saveTxtFl:
		{
			ret &= Temps.save2file();
			ret &= U00.save2file();
			ret &= U01.save2file();
			ret &= U10.save2file();
			ret &= U11.save2file();
			break;
		}
		case saveBinFl:
		{
			break;
		}
		}
	}
	return ret;
}

bool NodeV::saveFromDevice(bool saveAll = false, saveFlags saveFlag_ = saveFl,
	cl_command_queue * que_ = nullptr)
{
	bool ret = true;
	if (saveAll)
	{
		switch (saveFlag_)
		{
		case saveFl:
		{
		}
		case saveTxtFl:
		{
			ret &= Temps.save_txt_from_device("", que_);
			ret &= U00.save_txt_from_device("", que_);
			ret &= U01.save_txt_from_device("", que_);
			ret &= U10.save_txt_from_device("", que_);
			ret &= U11.save_txt_from_device("", que_);
			break;
		}
		case saveBinFl:
		{
			break;
		}
		}
	}
	return ret;
}


//////////////////////////////////////////////////////
///////////            BLinks            /////////////
//////////////////////////////////////////////////////
void BLinks::setStruct(legacyBLinks& struct_, const int i)
{
	vNvec(i) = struct_.vNvec;
	Tau(i) = struct_.Tau;
	blLen(i) = struct_.blLen;
	Node_loc(i) = struct_.Node_loc;
	P01ind(i) = struct_.P01ind;
	int_type(i) = struct_.int_type;
	if (glFlag)
		colorInds(i) = struct_.colorInd;
}

legacyBLinks BLinks::getStruct(const int i)
{
	legacyBLinks struct_;
	struct_.vNvec = vNvec(i);
	struct_.Tau = Tau(i);
	struct_.blLen = blLen(i);
	struct_.Node_loc = Node_loc(i);
	struct_.P01ind = P01ind(i);
	struct_.int_type = int_type(i);
	if (glFlag)
		struct_.colorInd = colorInds(i);
	else
		struct_.colorInd = 0;
	return struct_;
}


void BLinks::allocateArrays()
{
	vNvec.setName("BLinks_vN");
	vNvec.allocate(fullSize);
	vNvec.fill({ 0., 0. });

	Tau.setName("BLinks_Tau");
	Tau.allocate(fullSize);
	Tau.fill(0.);

	blLen.setName("BLinks_blLen");
	blLen.allocate(fullSize);
	blLen.fill(0.);

	Node_loc.setName("BLinks_nLoc");
	Node_loc.allocate(fullSize);
	Node_loc.fill(0);

	P01ind.setName("BLinks_vT");
	P01ind.allocate(fullSize);
	P01ind.fill({ 0, 0 });

	int_type.setName("BLinks_intType");
	int_type.allocate(fullSize);
	int_type.fill(0);
}

void BLinks::allocateBuffers(int bufSize = -1)
{
	if (bufSize == -1)
		bufSize = fullSize;

	bufFullSize = bufSize;

	vNvec.allocate_buffer_size(bufFullSize);
	Tau.allocate_buffer_size(bufFullSize);
	blLen.allocate_buffer_size(bufFullSize);
	Node_loc.allocate_buffer_size(bufFullSize);
	P01ind.allocate_buffer_size(bufFullSize);
	int_type.allocate_buffer_size(bufFullSize);
}

void BLinks::copyToDevice(int writeSize = -1, cl_command_queue * que_ = nullptr, cl_bool bFlag_ = CL_TRUE)
{
	if (writeSize == -1)
		writeSize = fullSize;

	vNvec.copy_to_buffer_size(writeSize, que_, bFlag_);
	Tau.copy_to_buffer_size(writeSize, que_, bFlag_);
	blLen.copy_to_buffer_size(writeSize, que_, bFlag_);
	Node_loc.copy_to_buffer_size(writeSize, que_, bFlag_);
	P01ind.copy_to_buffer_size(writeSize, que_, bFlag_);
	int_type.copy_to_buffer_size(writeSize, que_, bFlag_);
}

void BLinks::copyToHost(int readSize = -1, cl_command_queue * que_ = nullptr, cl_bool bFlag_ = CL_TRUE)
{
	if (readSize == -1)
		readSize = fullSize;

	vNvec.read_from_buffer_size(readSize, que_, bFlag_);
	Tau.read_from_buffer_size(readSize, que_, bFlag_);
	blLen.read_from_buffer_size(readSize, que_, bFlag_);
	Node_loc.read_from_buffer_size(readSize, que_, bFlag_);
	P01ind.read_from_buffer_size(readSize, que_, bFlag_);
	int_type.read_from_buffer_size(readSize, que_, bFlag_);
}

void BLinks::createColorInds()
{
	glFlag = true;

	colorInds.setName("BLinks_colorInds");
	colorInds.allocate(fullSize);
	colorInds.fill(0);

	colorInds.allocate_buffer_size(fullSize);
}

void BLinks::ini()
{
	//initialization functions specific to class
}

void BLinks::setBuffers(Kernel& ker, int& curind, arrName arrList[], int numArrs)
{
	for (int i = 0; i < numArrs; i++)
	{
		switch (arrList[i])
		{
		case vNArr:
		{
			ker.set_argument(curind, vNvec.get_buf_add());
			break;
		}
		case tauArr:
		{
			ker.set_argument(curind, Tau.get_buf_add());
			break;
		}
		case lenArr:
		{
			ker.set_argument(curind, blLen.get_buf_add());
			break;
		}
		case nodLocArr:
		{
			ker.set_argument(curind, Node_loc.get_buf_add());
			break;
		}
		case P01Arr:
		{
			ker.set_argument(curind, P01ind.get_buf_add());
			break;
		}
		case typeArr:
		{
			ker.set_argument(curind, int_type.get_buf_add());
			break;
		}
		case colorArr:
		{
			ker.set_argument(curind, colorInds.get_buf_add());
			break;
		}
		}
		curind++;
	}
}

bool BLinks::save(bool saveAll = false, saveFlags saveFlag_ = saveFl)
{
	bool ret = true;

	if (saveAll)
	{
		switch (saveFlag_)
		{
		case saveFl:
		{
		}
		case saveTxtFl:
		{
			ret &= vNvec.save2file();
			ret &= Tau.save2file();
			ret &= blLen.save2file();
			ret &= Node_loc.save2file();
			ret &= P01ind.save2file();
			ret &= int_type.save2file();
			if (glFlag)
				ret &= colorInds.save2file();
			break;
		}
		case saveBinFl:
		{
			break;
		}
		}
	}
	return ret;
}

bool BLinks::saveFromDevice(bool saveAll = false, saveFlags saveFlag_ = saveFl,
	cl_command_queue * que_ = nullptr)
{
	bool ret = true;
	if (saveAll)
	{
		switch (saveFlag_)
		{
		case saveFl:
		{
		}
		case saveTxtFl:
		{
			ret &= vNvec.save_txt_from_device("", que_);
			ret &= Tau.save_txt_from_device("", que_);
			ret &= blLen.save_txt_from_device("", que_);
			ret &= Node_loc.save_txt_from_device("", que_);
			ret &= P01ind.save_txt_from_device("", que_);
			ret &= int_type.save_txt_from_device("", que_);
			if (glFlag)
				ret &= colorInds.save_txt_from_device("",que_);
			break;
		}
		case saveBinFl:
		{
			break;
		}
		}
	}
	return ret;
}


//////////////////////////////////////////////////////
///////////            PParam            /////////////
//////////////////////////////////////////////////////