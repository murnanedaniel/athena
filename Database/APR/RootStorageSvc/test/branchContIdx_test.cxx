/*
  Copyright (C) 2002-2021 CERN for the benefit of the ATLAS collaboration
*/

/**
 * @file APR/RootStorageSvc/test/branchContIdx_test.cxx
 * @author Marcin Nowak
 * @brief Test BranchContainer in Branch Fill Mode with Indexing
 */


/*
  This test will write different number of objects (in this case DbStrings)
  (3 objects to Container A and 1 object to Container B to the same Database, in one transaction)
  and verify that the objects can be read back.
*/

#include "PersistentDataModel/Guid.h"
#include "PersistentDataModel/Token.h"

#include "StorageSvc/DbType.h"
#include "StorageSvc/Shape.h"
#include "StorageSvc/IStorageSvc.h"
#include "StorageSvc/IStorageExplorer.h"
#include "StorageSvc/FileDescriptor.h"
#include "StorageSvc/DbReflex.h"
#include "StorageSvc/DbString.h"

using namespace std;

const string            filename = "branchContIdx_testfile.root";    // test file name

// use normal containers
//const string            containerNameA = "ContainerA";               // container for objects A
//const string            containerNameB = "ContainerB";               // container for objects B

// use branch containes
const string            containerNameA = "Containers(A)";               // container for objects A
const string            containerNameB = "Containers(B)";               // container for objects B


//const pool::DbType      storageType = pool::ROOTTREE_StorageType;
const pool::DbType      storageType = pool::ROOTTREEINDEX_StorageType;

pool::DbString       myStringA1("This is my string A1");
pool::DbString       myStringA2("This is my string A2");
pool::DbString       myStringA3("This is my string A3");
pool::DbString       myStringA4("This is my string A4");
pool::DbString       myStringA5("This is my string A5");
pool::DbString       myStringA6("This is my string A6");

pool::DbString       myStringB1("This is my string B1");
pool::DbString       myStringB2("This is my string B2");

// string tokens to written objects, for use when reading back
string refA1, refA2, refA3, refA4, refA5, refA6, refB1, refB2;


int main() {
   cout << "Starting branchContIdx_test" << endl;

   pool::IStorageSvc* storSvc = pool::createStorageSvc("StorageSvc");
   if ( ! storSvc ) {
      throw std::runtime_error( "Could not create a StorageSvc object" );
   }
   storSvc->addRef();
   cout << "Start WRITE session" << endl;
   pool::Session* sessionHandle = 0;
   if ( ! ( storSvc->startSession( pool::RECREATE, storageType.type(), sessionHandle ).isSuccess() ) ) {
      throw std::runtime_error( "Could not start a session." );
   }

   cout << "Session connect" << endl;
   pool::FileDescriptor fd( filename, filename );
   if ( ! ( storSvc->connect( sessionHandle, pool::RECREATE, fd ).isSuccess() ) ) {
      throw std::runtime_error( "Could not start a connection." );
   }
   pool::DatabaseConnection* connection = fd.dbc();

   // Retrieve DbString dictionary
   RootType class_String ( "pool::DbString" );
   if ( ! class_String ) {
      throw std::runtime_error( "Could not retrieve the dictionary for class DbString" );
   }
   // Create shape for DbString
   Guid guid = pool::DbReflex::guid(class_String);
   const pool::Shape* shape = 0;
   if ( storSvc->getShape( fd, guid, shape ) == pool::IStorageSvc::SHAPE_NOT_AVAILIBLE ) {
      cout << "need to create a Shape for DbString" << endl;
      storSvc->createShape( fd, containerNameA, guid, shape );
   }
   if( ! shape ) {
      throw std::runtime_error( "Could not create a persistent shape." );
   }

   // Commit here to test empty commits
   if( ! ( storSvc->endTransaction( connection, pool::Transaction::TRANSACT_COMMIT ).isSuccess() ) ) {
      throw std::runtime_error( "Empty commit FAILED" );
   }

   // helper lambda to write a string object into a container
   auto writeStr = [&](const string& contName, pool::DbString& str) {
      Token* token{nullptr};
      if( !storSvc->allocate( fd, contName, storageType.type(), &str, shape, token ).isSuccess() ) {
         throw std::runtime_error( "Could not write an object" );
      }
      return token;
   };
   
   refA1 = writeStr(containerNameA, myStringA1)->toString();
   refA2 = writeStr(containerNameA, myStringA2)->toString();
   refA3 = writeStr(containerNameA, myStringA3)->toString();
   refB1 = writeStr(containerNameB, myStringB1)->toString();

   // Commit
   if( ! ( storSvc->endTransaction( connection, pool::Transaction::TRANSACT_COMMIT ).isSuccess() ) ) {
      throw std::runtime_error( "Commit FAILED" );
   }

   refA4 = writeStr(containerNameA, myStringA4)->toString();
   refB2 = writeStr(containerNameB, myStringB2)->toString();
   refA5 = writeStr(containerNameA, myStringA5)->toString();
   refA6 = writeStr(containerNameA, myStringA6)->toString();

   // Commit again
   if( ! ( storSvc->endTransaction( connection, pool::Transaction::TRANSACT_COMMIT ).isSuccess() ) ) {
      throw std::runtime_error( "Commit FAILED" );
   }

   // Close session
   if( !storSvc->disconnect( fd ).isSuccess() ) {
      throw std::runtime_error( "Could not disconnect." );
   }
   if( !storSvc->endSession( sessionHandle ).isSuccess() ) {
      throw std::runtime_error( "Could not end correctly the session." );
   }

   // not realeasing if reading back in the same function
   // storSvc->release();

   cout << "Sample Tokens to written objects:" << endl;
   cout << "Token for object A1 = " << refA1 << endl;
   cout << "Token for object A2 = " << refA2 << endl;
   cout << "Token for object A3 = " << refA3 << endl;
   cout << "Token for object A4 = " << refA4 << endl;
   cout << "Token for object A5 = " << refA5 << endl;
   cout << "Token for object A6 = " << refA6 << endl;
   cout << "Token for object B1 = " << refB1 << endl;
   cout << "Token for object B2 = " << refB2 << endl;


   // ===============    READ back
   cout << endl << "Starting READ" << endl;
   sessionHandle = nullptr;
   if( !storSvc->startSession( pool::READ, storageType.type(), sessionHandle ).isSuccess() ) {
      throw std::runtime_error( "Could not start the read session." );
   }
   if( storSvc->connect( sessionHandle, pool::READ, fd ) != pool::DbStatus::Success ) {
      throw std::runtime_error( "Could not start a read connection." );
   }
   // get shape again
   shape = nullptr;
   if ( storSvc->getShape( fd, guid, shape ) == pool::IStorageSvc::SHAPE_NOT_AVAILIBLE ) {
      cout << "need to create a Shape for DbString" << endl;
      storSvc->createShape( fd, containerNameA, guid, shape );
   }
   if( ! shape ) {
      throw std::runtime_error( "Could not create a persistent shape." );
   }

   pool::DbString readString;
   auto readfun = [&](const string& ref) {
      Token token;
      token.fromString( ref );
      pool::DbString *str_p = &readString;
      if( !storSvc->read(fd, token, shape, (void**)&str_p).isSuccess() ) {
         throw std::runtime_error( "Read FAILED" );
      }
   };
   
   readfun(refA1);   cout << "read back A1: " << readString << endl;
   readfun(refB1);   cout << "read back B1: " << readString << endl;
   readfun(refA3);   cout << "read back A3: " << readString << endl;
   readfun(refA6);   cout << "read back A6: " << readString << endl;
   readfun(refA4);   cout << "read back A4: " << readString << endl;
   readfun(refB2);   cout << "read back B2: " << readString << endl;
   readfun(refA5);   cout << "read back A5: " << readString << endl;
   readfun(refA2);   cout << "read back A2: " << readString << endl;
   
   // Close session
   if( !storSvc->disconnect( fd ).isSuccess() ) {
      throw std::runtime_error( "Could not disconnect." );
   }
   if( !storSvc->endSession( sessionHandle ).isSuccess() ) {
      throw std::runtime_error( "Could not end correctly the session." );
   }

   storSvc->release();

   return 0;
}

